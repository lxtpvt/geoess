#' Fit a Multivariate Normal Mixture by EM
#'
#' Fits a mixture of multivariate Normals to parameter draws.
#' EM structure is inspired by standard mixture fitting workflows (e.g., RBesT),
#' but implemented independently for this package.
#'
#' @param draws Numeric matrix of draws (iterations x parameters).
#' @param family Character. "mvnorm" or "norm".
#' @param K Integer. Number of mixture components.
#' @param init Character. Initialization method, "kmeans" or "random".
#' @param max_iter Integer. Maximum EM iterations.
#' @param tol Numeric. Convergence tolerance on log-likelihood.
#' @param min_weight Numeric. Minimum component weight.
#' @param ridge Numeric. Diagonal ridge added to covariances.
#'
#' @return A list of class \code{geoess_mix} with mixture parameters and fit diagnostics.
#' @export
geoess_mixfit <- function(draws, family = c("mvnorm", "norm"), K,
                          init = c("kmeans", "random"), max_iter = 500,
                          tol = 1e-8, min_weight = 1e-8, ridge = 1e-8) {
  family <- match.arg(family)
  init <- match.arg(init)
  if (missing(K) || length(K) != 1L || K <= 0) stop("K must be a positive integer")

  X <- as.matrix(draws)
  if (!is.numeric(X)) stop("draws must be numeric")
  if (any(!is.finite(X))) stop("draws contain non-finite values")

  n <- nrow(X)
  d <- ncol(X)
  if (n < K) stop("Number of draws must be >= K")

  if (family == "norm" && d != 1L) {
    stop("family = 'norm' requires univariate draws")
  }

  if (K == 1L) {
    mu <- colMeans(X)
    if (n < 2L) {
      Sigma <- diag(ridge, d)
    } else {
      Sigma <- stats::cov(X)
      if (d == 1L) Sigma <- matrix(as.numeric(Sigma), nrow = 1L)
      Sigma <- Sigma + diag(ridge, d)
    }
    if (any(!is.finite(Sigma))) Sigma <- diag(ridge, d)

    logdens <- mvtnorm::dmvnorm(X, mean = mu, sigma = Sigma, log = TRUE)
    logLik <- sum(logdens)
    n_params <- (K - 1) + K * d + K * d * (d + 1) / 2
    out <- list(
      family = "mvnorm",
      K = K,
      weights = 1,
      mu = list(mu),
      Sigma = list(Sigma),
      prec = list(.geoess_solve(Sigma)),
      chol = list(.geoess_try_chol(Sigma)$R),
      logLik = logLik,
      AIC = -2 * logLik + 2 * n_params,
      BIC = -2 * logLik + log(n) * n_params,
      iter = 1L,
      converged = TRUE
    )
    class(out) <- "geoess_mix"
    return(out)
  }

  if (init == "kmeans") {
    km <- stats::kmeans(X, centers = K, nstart = 5)
    mu <- lapply(seq_len(K), function(k) km$centers[k, ])
    weights <- tabulate(km$cluster, nbins = K) / n
  } else {
    idx <- sample(seq_len(n), K)
    mu <- lapply(idx, function(i) X[i, , drop = FALSE])
    mu <- lapply(mu, function(v) as.numeric(v))
    weights <- rep(1 / K, K)
  }

  Sigma <- vector("list", K)
  for (k in seq_len(K)) {
    if (init == "kmeans") {
      cluster_idx <- which(km$cluster == k)
      if (length(cluster_idx) > 1L) {
        Sigma[[k]] <- stats::cov(X[cluster_idx, , drop = FALSE])
      } else {
        Sigma[[k]] <- stats::cov(X)
      }
    } else {
      Sigma[[k]] <- stats::cov(X)
    }
    if (d == 1L) Sigma[[k]] <- matrix(as.numeric(Sigma[[k]]), nrow = 1L)
    Sigma[[k]] <- Sigma[[k]] + diag(ridge, d)
    if (any(!is.finite(Sigma[[k]]))) Sigma[[k]] <- diag(ridge, d)
  }

  logLik <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    log_comp <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in seq_len(K)) {
      log_comp[, k] <- log(weights[k]) + mvtnorm::dmvnorm(
        X, mean = mu[[k]], sigma = Sigma[[k]], log = TRUE
      )
    }

    log_mix <- apply(log_comp, 1L, .geoess_log_sum_exp)
    resp <- exp(log_comp - log_mix)
    logLik_new <- sum(log_mix)

    if (is.finite(logLik) && abs(logLik_new - logLik) < tol * (1 + abs(logLik))) {
      converged <- TRUE
      logLik <- logLik_new
      break
    }

    logLik <- logLik_new

    Nk <- colSums(resp)
    weights <- pmax(Nk / n, min_weight)
    weights <- weights / sum(weights)

    for (k in seq_len(K)) {
      if (Nk[k] <= 1e-8) {
        mu[[k]] <- mu[[k]]
        Sigma[[k]] <- Sigma[[k]] + diag(ridge, d)
        next
      }

      mu[[k]] <- as.numeric(colSums(resp[, k] * X) / Nk[k])
      centered <- sweep(X, 2, mu[[k]], FUN = "-")
      weighted <- centered * sqrt(resp[, k])
      Sigma[[k]] <- crossprod(weighted) / Nk[k] + diag(ridge, d)
    }
  }

  n_params <- (K - 1) + K * d + K * d * (d + 1) / 2
  out <- list(
    family = "mvnorm",
    K = K,
    weights = weights,
    mu = mu,
    Sigma = Sigma,
    prec = lapply(Sigma, .geoess_solve),
    chol = lapply(Sigma, function(S) .geoess_try_chol(S)$R),
    logLik = logLik,
    AIC = -2 * logLik + 2 * n_params,
    BIC = -2 * logLik + log(n) * n_params,
    iter = iter,
    converged = converged
  )
  class(out) <- "geoess_mix"
  out
}

#' Automated Mixture Selection
#'
#' Fits mixtures across K and selects by AIC/BIC.
#'
#' @param draws Numeric matrix of draws.
#' @param family Character. "mvnorm" or "norm".
#' @param K_grid Integer vector of candidate component counts.
#' @param criterion Character. "BIC" or "AIC".
#' @param ... Passed to \code{geoess_mixfit}.
#'
#' @return A fitted \code{geoess_mix} object with selection metadata.
#' @export
geoess_automixfit <- function(draws, family = c("mvnorm", "norm"),
                              K_grid = 1:6, criterion = c("BIC", "AIC"), ...) {
  family <- match.arg(family)
  criterion <- match.arg(criterion)
  if (!is.numeric(K_grid) || length(K_grid) == 0L) {
    stop("K_grid must be a non-empty numeric vector")
  }

  fits <- list()
  crit_vals <- numeric(length(K_grid))

  for (i in seq_along(K_grid)) {
    K <- K_grid[i]
    fit <- tryCatch(
      geoess_mixfit(draws, family = family, K = K, ...),
      error = function(e) NULL
    )
    fits[[i]] <- fit
    if (is.null(fit)) {
      crit_vals[i] <- Inf
    } else {
      crit_vals[i] <- fit[[criterion]]
    }
  }

  best_idx <- which.min(crit_vals)
  if (is.infinite(crit_vals[best_idx])) stop("All mixture fits failed")

  best_fit <- fits[[best_idx]]
  best_fit$criterion <- criterion
  best_fit$K_selected <- K_grid[best_idx]
  best_fit$criterion_values <- crit_vals
  best_fit
}
