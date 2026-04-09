#' Schur Complement Helper (Internal)
#'
#' Computes the Schur complement for a block of a square matrix.
#'
#' @param M Numeric matrix.
#' @param idx_a Integer indices for the target block.
#' @param idx_b Integer indices for the conditioning block. If NULL, uses
#'   the complement of \code{idx_a}.
#' @param symmetrize Logical; if TRUE, symmetrize the result.
#' @param solver Character. "chol" (preferred) or "solve".
#'
#' @return Numeric matrix for the conditional block.
#' @keywords internal
.geoess_schur <- function(M, idx_a, idx_b = NULL,
                          symmetrize = TRUE,
                          solver = c("chol", "solve")) {
  solver <- match.arg(solver)
  M <- .geoess_as_matrix(M, "M")
  .geoess_check_square(M, "M")

  d <- nrow(M)
  if (is.null(idx_a) || length(idx_a) == 0L) stop("idx_a must be non-empty")
  idx_a <- as.integer(idx_a)
  if (any(!is.finite(idx_a))) stop("idx_a must be finite integers")
  if (any(idx_a < 1L | idx_a > d)) stop("idx_a out of bounds")
  if (any(duplicated(idx_a))) stop("idx_a contains duplicates")

  if (is.null(idx_b)) {
    idx_b <- setdiff(seq_len(d), idx_a)
  } else {
    idx_b <- as.integer(idx_b)
    if (length(idx_b) == 0L) {
      return(M[idx_a, idx_a, drop = FALSE])
    }
    if (any(!is.finite(idx_b))) stop("idx_b must be finite integers")
    if (any(idx_b < 1L | idx_b > d)) stop("idx_b out of bounds")
    if (any(duplicated(idx_b))) stop("idx_b contains duplicates")
    if (length(intersect(idx_a, idx_b)) > 0L) {
      stop("idx_a and idx_b must be disjoint")
    }
  }

  if (length(idx_b) == 0L) {
    out <- M[idx_a, idx_a, drop = FALSE]
    return(if (symmetrize) .geoess_symmetrize(out) else out)
  }

  M_aa <- M[idx_a, idx_a, drop = FALSE]
  M_ab <- M[idx_a, idx_b, drop = FALSE]
  M_bb <- M[idx_b, idx_b, drop = FALSE]
  M_ba <- M[idx_b, idx_a, drop = FALSE]

  tmp <- NULL
  if (solver == "chol") {
    chol_fit <- .geoess_try_chol(as.matrix(M_bb))
    if (!is.null(chol_fit$R)) {
      R <- chol_fit$R
      tmp <- backsolve(R, forwardsolve(t(R), M_ba))
    } else {
      warning("Block matrix is not positive definite; using solve fallback")
    }
  }
  if (is.null(tmp)) {
    tmp <- .geoess_solve(M_bb, M_ba)
  }

  out <- M_aa - M_ab %*% tmp
  if (symmetrize) out <- .geoess_symmetrize(out)
  out
}

#' Conditional ESS via Schur Complement
#'
#' Computes conditional ESS for a component (or block) of parameters when
#' parameters are correlated, using Schur complements on \code{I1} and \code{Ipi}.
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param idx_a Integer indices for the target component/block.
#' @param idx_b Integer indices for the conditioning component/block.
#' @param method Character. "trace" or "scalar".
#' @param symmetrize Logical; if TRUE, symmetrize intermediate results.
#' @param solver Character. "chol" or "solve".
#'
#' @return A list with \code{I1_cond}, \code{Ipi_cond}, and \code{ess_cond}.
#' @export
geoess_conditional <- function(I1, Ipi, idx_a, idx_b,
                               method = c("trace", "scalar"),
                               symmetrize = TRUE,
                               solver = c("chol", "solve")) {
  method <- match.arg(method)
  solver <- match.arg(solver)

  I1 <- .geoess_as_matrix(I1, "I1")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I1, "I1")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I1, Ipi, "I1", "Ipi")

  I1_cond <- .geoess_schur(I1, idx_a, idx_b, symmetrize = symmetrize, solver = solver)
  Ipi_cond <- .geoess_schur(Ipi, idx_a, idx_b, symmetrize = symmetrize, solver = solver)

  if (method == "scalar") {
    if (nrow(I1_cond) != 1L || ncol(I1_cond) != 1L) {
      stop("method = 'scalar' is only valid for a scalar component")
    }
    ess_cond <- as.numeric(Ipi_cond) / as.numeric(I1_cond)
  } else {
    ess_cond <- geoess_trace(I1_cond, Ipi_cond, symmetrize = symmetrize)
  }

  list(I1_cond = I1_cond, Ipi_cond = Ipi_cond, ess_cond = ess_cond)
}

#' Conditional ESS at a Reference Point
#'
#' @param I1_fn Fisher information function or matrix.
#' @param Ipi_fn Prior curvature function, matrix, or \code{geoess_mix} object.
#' @param theta_ref Numeric vector. Reference point.
#' @param idx_a Integer indices for the target component/block.
#' @param idx_b Integer indices for the conditioning component/block.
#' @param ... Passed to \code{I1_fn} or \code{Ipi_fn}.
#'
#' @return A list with conditional matrices and ESS.
#' @export
geoess_conditional_local <- function(I1_fn, Ipi_fn, theta_ref,
                                     idx_a, idx_b, ...) {
  theta_ref <- as.numeric(theta_ref)

  I1_val <- if (is.function(I1_fn)) I1_fn(theta_ref, ...) else I1_fn
  Ipi_val <- if (inherits(Ipi_fn, "geoess_mix")) {
    geoess_Ipi_from_mix(Ipi_fn, theta_ref)
  } else if (is.function(Ipi_fn)) {
    Ipi_fn(theta_ref, ...)
  } else {
    Ipi_fn
  }

  geoess_conditional(I1 = I1_val, Ipi = Ipi_val, idx_a = idx_a, idx_b = idx_b, ...)
}

#' Global Conditional ESS under a Weighting Distribution
#'
#' Computes the expectation of conditional ESS under a weighting distribution.
#'
#' @param Q Distribution for \eqn{\theta}. See \code{geoess_global}.
#' @param I1 Fisher information per observation.
#' @param Ipi Prior curvature.
#' @param idx_a Integer indices for the target component/block.
#' @param idx_b Integer indices for the conditioning component/block.
#' @param nsim Integer. Number of draws if sampling is needed.
#' @param seed Integer. Seed for reproducibility.
#' @param weights Optional weights aligned with draws.
#' @param psd Character. "none", "clip", or "nearPD" for Ipi adjustment.
#' @param return_samples Logical; if TRUE return per-draw ESS values.
#' @param ... Passed to \code{I1} or \code{Ipi} if they are functions.
#'
#' @return A list with estimated ESS and uncertainty summaries.
#' @export
geoess_conditional_global <- function(Q,
                                      I1,
                                      Ipi,
                                      idx_a,
                                      idx_b,
                                      nsim = 2000,
                                      seed = 1,
                                      weights = NULL,
                                      psd = c("none", "clip", "nearPD"),
                                      return_samples = FALSE,
                                      ...) {
  psd <- match.arg(psd)
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim <= 0) {
    stop("nsim must be a positive integer")
  }

  q <- .geoess_Q_draws(Q, nsim = nsim, seed = seed)
  draws <- q$draws
  n <- nrow(draws)
  d <- ncol(draws)
  if (n == 0) stop("No draws available for Q")
  if (is.null(d) || d < 1L) stop("Draws must have at least one column")

  w_all <- if (!is.null(q$weights)) q$weights else NULL
  if (!is.null(weights)) {
    if (is.null(w_all)) {
      w_all <- weights
    } else {
      if (length(weights) != length(w_all)) stop("weights length mismatch")
      w_all <- w_all * weights
    }
  }
  w_norm <- .geoess_normalize_weights(w_all, n)
  n_eff <- .geoess_is_ess(w_norm)

  I1_get <- .geoess_I1_getter(I1, n, d, ...)
  Ipi_get <- .geoess_Ipi_getter(Ipi, n, d, psd, ...)
  I1_scalar <- is.numeric(I1) && length(I1) == 1L

  ess_vec <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    theta <- draws[i, ]
    I1_i <- I1_get(i, theta)
    if (I1_scalar && is.matrix(I1_i) && (nrow(I1_i) == 0 || ncol(I1_i) == 0)) {
      d_loc <- length(theta)
      if (d_loc < 1L) stop("I1 dimension is invalid")
      I1_i <- .geoess_diag_scalar(as.numeric(I1), d_loc)
    }
    Ipi_i <- Ipi_get(i, theta)
    if (any(!is.finite(as.numeric(I1_i))) || any(!is.finite(as.numeric(Ipi_i)))) {
      next
    }

    cond_res <- geoess_conditional(I1_i, Ipi_i, idx_a = idx_a, idx_b = idx_b)
    ess_vec[i] <- as.numeric(cond_res$ess_cond)
  }

  valid <- is.finite(ess_vec)
  if (!all(valid)) {
    warning("Dropping draws with non-finite conditional ESS values")
    ess_vec <- ess_vec[valid]
    w_norm <- w_norm[valid]
    w_norm <- w_norm / sum(w_norm)
    n_eff <- .geoess_is_ess(w_norm)
  }
  if (length(ess_vec) == 0) stop("No valid conditional ESS samples")

  estimate <- .geoess_weighted_mean(ess_vec, w_norm)
  mcse <- sqrt(.geoess_weighted_var(ess_vec, w_norm)) / sqrt(n_eff)
  quant_probs <- c(0.05, 0.5, 0.95)
  quant <- .geoess_weighted_quantile(ess_vec, w_norm, probs = quant_probs)
  names(quant) <- paste0("q", quant_probs)

  res <- list(
    type = if (!is.null(q$type)) q$type else "custom",
    method = "conditional",
    estimate = estimate,
    mcse = mcse,
    n = length(ess_vec),
    quantiles = quant,
    weights_summary = if (!is.null(w_all)) .geoess_weight_summary(w_norm) else NULL,
    ess_samples = if (return_samples) ess_vec else NULL,
    call = match.call()
  )
  class(res) <- "geoess_global_result"
  res
}

#' ESS Spectrum with Eigenvectors
#'
#' Returns eigenvalues and eigenvectors of the symmetrized operator
#' \eqn{M = I_1^{-1} I_\pi}.
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param symmetrize Logical; if TRUE, symmetrize \eqn{M}.
#'
#' @return A list with \code{values}, \code{vectors}, and \code{M}.
#' @export
geoess_spectrum <- function(I1, Ipi, symmetrize = TRUE) {
  I1 <- .geoess_as_matrix(I1, "I1")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I1, "I1")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I1, Ipi, "I1", "Ipi")

  M <- .geoess_solve(I1, Ipi)
  if (symmetrize) M <- .geoess_symmetrize(M)

  eig <- eigen(M, symmetric = TRUE)
  list(values = eig$values, vectors = eig$vectors, M = M)
}
