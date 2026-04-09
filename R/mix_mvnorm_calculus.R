#' MVN Mixture Log-Density
#'
#' Computes log q(theta) for a finite MVN mixture.
#'
#' @param mix A \code{geoess_mix} object with MVN components.
#' @param theta Numeric vector evaluation point.
#' @param use_cache Logical; if TRUE, uses cached Cholesky factors if present.
#'
#' @return Numeric scalar log density.
#' @export
geoess_mix_logdens <- function(mix, theta, use_cache = TRUE) {
  if (!inherits(mix, "geoess_mix")) stop("mix must be a geoess_mix object")
  if (is.null(mix$family) || mix$family != "mvnorm") {
    stop("geoess_mix_logdens currently supports mix$family == 'mvnorm'")
  }

  theta <- as.numeric(theta)
  if (!all(is.finite(theta))) stop("theta contains non-finite values")

  w <- as.numeric(mix$weights)
  if (any(!is.finite(w)) || any(w < 0)) stop("mix$weights must be finite and nonnegative")
  w <- w / sum(w)

  K <- length(w)
  log_comp <- numeric(K)
  for (k in seq_len(K)) {
    mu_k <- as.numeric(mix$mu[[k]])
    Sig_k <- mix$Sigma[[k]]
    chol_k <- NULL
    if (use_cache && !is.null(mix$chol) && length(mix$chol) >= k) {
      chol_k <- mix$chol[[k]]
    }
    log_comp[k] <- log(w[k]) + .geoess_log_mvn(theta, mu_k, Sig_k, chol = chol_k)
  }
  .geoess_logsumexp(log_comp)
}

#' MVN Mixture Gradient of Log-Density
#'
#' @param mix A \code{geoess_mix} object with MVN components.
#' @param theta Numeric vector evaluation point.
#' @param use_cache Logical; if TRUE, uses cached precision/Cholesky if present.
#'
#' @return Numeric vector gradient.
#' @export
geoess_mix_gradlog <- function(mix, theta, use_cache = TRUE) {
  if (!inherits(mix, "geoess_mix")) stop("mix must be a geoess_mix object")
  if (is.null(mix$family) || mix$family != "mvnorm") {
    stop("geoess_mix_gradlog currently supports mix$family == 'mvnorm'")
  }

  theta <- as.numeric(theta)
  if (!all(is.finite(theta))) stop("theta contains non-finite values")

  w <- as.numeric(mix$weights)
  if (any(!is.finite(w)) || any(w < 0)) stop("mix$weights must be finite and nonnegative")
  w <- w / sum(w)

  K <- length(w)
  d <- length(theta)
  log_comp <- numeric(K)
  for (k in seq_len(K)) {
    mu_k <- as.numeric(mix$mu[[k]])
    Sig_k <- mix$Sigma[[k]]
    chol_k <- NULL
    if (use_cache && !is.null(mix$chol) && length(mix$chol) >= k) {
      chol_k <- mix$chol[[k]]
    }
    log_comp[k] <- log(w[k]) + .geoess_log_mvn(theta, mu_k, Sig_k, chol = chol_k)
  }
  resp <- .geoess_softmax_log(log_comp)

  g <- numeric(d)
  for (k in seq_len(K)) {
    mu_k <- as.numeric(mix$mu[[k]])
    prec_k <- NULL
    if (use_cache && !is.null(mix$prec) && length(mix$prec) >= k && !is.null(mix$prec[[k]])) {
      prec_k <- mix$prec[[k]]
    } else if (use_cache && !is.null(mix$chol) && length(mix$chol) >= k && !is.null(mix$chol[[k]])) {
      prec_k <- .geoess_prec_from_chol(mix$chol[[k]])
    } else {
      prec_k <- .geoess_solve(mix$Sigma[[k]])
    }
    gk <- -as.numeric(prec_k %*% (theta - mu_k))
    g <- g + resp[k] * gk
  }
  g
}

#' MVN Mixture Hessian of Log-Density
#'
#' @param mix A \code{geoess_mix} object with MVN components.
#' @param theta Numeric vector evaluation point.
#' @param use_cache Logical; if TRUE, uses cached precision/Cholesky if present.
#' @param symmetrize Logical; if TRUE, symmetrize the result.
#'
#' @return Numeric matrix Hessian of log q(theta).
#' @export
geoess_mix_hesslog <- function(mix, theta, use_cache = TRUE, symmetrize = TRUE) {
  if (!inherits(mix, "geoess_mix")) stop("mix must be a geoess_mix object")
  if (is.null(mix$family) || mix$family != "mvnorm") {
    stop("geoess_mix_hesslog currently supports mix$family == 'mvnorm'")
  }

  theta <- as.numeric(theta)
  if (!all(is.finite(theta))) stop("theta contains non-finite values")

  w <- as.numeric(mix$weights)
  if (any(!is.finite(w)) || any(w < 0)) stop("mix$weights must be finite and nonnegative")
  w <- w / sum(w)

  K <- length(w)
  d <- length(theta)
  log_comp <- numeric(K)
  for (k in seq_len(K)) {
    mu_k <- as.numeric(mix$mu[[k]])
    Sig_k <- mix$Sigma[[k]]
    chol_k <- NULL
    if (use_cache && !is.null(mix$chol) && length(mix$chol) >= k) {
      chol_k <- mix$chol[[k]]
    }
    log_comp[k] <- log(w[k]) + .geoess_log_mvn(theta, mu_k, Sig_k, chol = chol_k)
  }
  resp <- .geoess_softmax_log(log_comp)

  g <- numeric(d)
  sum_H <- matrix(0, d, d)

  for (k in seq_len(K)) {
    mu_k <- as.numeric(mix$mu[[k]])
    prec_k <- NULL
    if (use_cache && !is.null(mix$prec) && length(mix$prec) >= k && !is.null(mix$prec[[k]])) {
      prec_k <- mix$prec[[k]]
    } else if (use_cache && !is.null(mix$chol) && length(mix$chol) >= k && !is.null(mix$chol[[k]])) {
      prec_k <- .geoess_prec_from_chol(mix$chol[[k]])
    } else {
      prec_k <- .geoess_solve(mix$Sigma[[k]])
    }

    diff <- theta - mu_k
    gk <- -as.numeric(prec_k %*% diff)
    Hk <- -prec_k
    g <- g + resp[k] * gk
    sum_H <- sum_H + resp[k] * (Hk + tcrossprod(gk))
  }

  H <- sum_H - tcrossprod(g)
  if (symmetrize) H <- .geoess_symmetrize(H)
  H
}

#' Prior Curvature from an MVN Mixture (Fisher-geometry ready)
#'
#' Computes the prior curvature matrix Ipi(theta) = - Hessian_theta log q(theta),
#' where q is a finite mixture of multivariate normal components.
#'
#' For q(theta) = sum_k w_k N(theta | mu_k, Sigma_k), define responsibilities
#' r_k(theta) proportional to w_k N(theta | mu_k, Sigma_k).
#' Let g_k = grad log N_k and H_k = Hess log N_k. Then:
#'   grad log q = sum_k r_k g_k
#'   Hess log q = sum_k r_k (H_k + g_k g_k^T) - (grad log q)(grad log q)^T
#' and Ipi = - Hess log q.
#'
#' @param mix A \code{geoess_mix} object with MVN components. Must contain:
#'   \code{weights} numeric length-K; \code{mu} list length-K of length-d vectors;
#'   \code{Sigma} list length-K of dxd cov matrices.
#'   Optional caches: \code{prec} list of dxd precision matrices or \code{chol} list of chol(Sigma).
#' @param theta Numeric vector (length d) where curvature is evaluated.
#' @param use_cache Logical; if TRUE, uses \code{mix$prec} or \code{mix$chol} if present.
#' @param symmetrize Logical; if TRUE, symmetrizes result for numerical stability.
#' @param psd Character. "none", "clip", or "nearPD" to enforce PSD.
#'
#' @return Numeric matrix (d x d), the prior curvature Ipi(theta).
#' @export
geoess_Ipi_from_mix <- function(mix, theta, use_cache = TRUE, symmetrize = TRUE,
                                psd = c("none", "clip", "nearPD")) {
  psd <- match.arg(psd)
  if (!inherits(mix, "geoess_mix")) stop("mix must be a geoess_mix object")
  if (is.null(mix$family) || mix$family != "mvnorm") {
    stop("geoess_Ipi_from_mix currently supports mix$family == 'mvnorm'")
  }

  theta <- as.numeric(theta)
  if (!all(is.finite(theta))) stop("theta contains non-finite values")

  w <- as.numeric(mix$weights)
  if (any(!is.finite(w)) || any(w < 0)) stop("mix$weights must be finite and nonnegative")
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("mix$weights must sum to a positive value")
  w <- w / sw

  K <- length(w)
  if (!is.list(mix$mu) || !is.list(mix$Sigma)) stop("mix$mu and mix$Sigma must be lists")
  if (length(mix$mu) != K || length(mix$Sigma) != K) {
    stop("mix$mu and mix$Sigma must have same length as mix$weights")
  }

  d <- length(theta)
  if (length(mix$mu[[1]]) != d) stop("theta dimension does not match mixture parameters")

  H <- geoess_mix_hesslog(mix, theta, use_cache = use_cache, symmetrize = FALSE)

  # Prior curvature: Ipi = - Hess log q
  Ipi <- -H
  if (symmetrize) Ipi <- .geoess_symmetrize(Ipi)

  if (psd == "clip") {
    ev <- eigen(Ipi, symmetric = TRUE)
    ev$values[ev$values < 0] <- 0
    Ipi <- ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
  } else if (psd == "nearPD") {
    Ipi <- as.matrix(Matrix::nearPD(Ipi, ensureSymmetry = TRUE)$mat)
  }
  Ipi
}
