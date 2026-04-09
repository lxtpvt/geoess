#' Beta Mixture Distributions
#'
#' Create and evaluate finite mixtures of Beta distributions.
#'
#' @param ... Components as numeric vectors of length 3 or a 3xK matrix.
#'   Each component is \code{c(weight, p1, p2)} where \code{p1, p2} are
#'   parameterized by \code{param}.
#' @param param Character. "ab" for shape1/shape2, "ms" for mean/sd,
#'   or "mn" for mean/effective sample size.
#'
#' @return A matrix of class \code{c("betaMix", "mix")} with rows
#'   \code{w}, \code{a}, \code{b}.
#' @export
mixbeta <- function(..., param = c("ab", "ms", "mn")) {
  mix <- .geoess_mixdist3(...)
  param <- match.arg(param)

  if (param == "ab") {
    ab <- t(mix[2:3, , drop = FALSE])
  } else if (param == "ms") {
    ab <- .geoess_ms2beta(mix[2, ], mix[3, ], drop = FALSE)
  } else {
    ab <- .geoess_mn2beta(mix[2, ], mix[3, ], drop = FALSE)
  }

  mix[2, ] <- ab[, 1]
  mix[3, ] <- ab[, 2]
  rownames(mix) <- c("w", "a", "b")

  if (any(mix["a", ] < 0) || any(mix["b", ] < 0)) {
    stop("beta shapes must be nonnegative")
  }

  class(mix) <- c("betaMix", "mix")
  attr(mix, "likelihood") <- "binomial"
  mix
}

#' Beta Prior Curvature
#'
#' Computes the local curvature of a Beta(a, b) prior:
#' \deqn{I_\pi(p) = (a-1)/p^2 + (b-1)/(1-p)^2}.
#' This is the curvature on the probability scale \eqn{p \in (0, 1)} and is
#' not the canonical-logit curvature used by
#' \code{\link{geoess_beta_binomial_canonical}}.
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @param a Numeric shape1 parameter (> 0).
#' @param b Numeric shape2 parameter (> 0).
#' @param eps Numeric. Clamp p to \code{[eps, 1-eps]} for stability.
#'
#' @return A 1x1 matrix (if length(p) == 1) or a diagonal matrix with the
#'   curvature values on the diagonal.
#' @export
Ipi_beta <- function(p, a, b, eps = 1e-6) {
  p <- as.numeric(p)
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L) {
    stop("a and b must be numeric scalars")
  }
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("a and b must be positive and finite")
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("eps must be in (0, 0.5)")
  }
  p <- pmin(pmax(p, eps), 1 - eps)
  vals <- (a - 1) / p^2 + (b - 1) / (1 - p)^2
  if (length(vals) == 1L) {
    return(matrix(as.numeric(vals), nrow = 1L, ncol = 1L))
  }
  diag(as.numeric(vals), nrow = length(vals))
}

#' @export
dmix <- function(mix, x, log = FALSE) {
  UseMethod("dmix")
}

#' @method dmix betaMix
#' @export
dmix.betaMix <- function(mix, x, log = FALSE) {
  .geoess_check_betaMix(mix)
  x <- as.numeric(x)
  w <- mix["w", ]
  a <- mix["a", ]
  b <- mix["b", ]
  K <- length(w)
  Nx <- length(x)

  log_comp <- matrix(NA_real_, nrow = Nx, ncol = K)
  for (k in seq_len(K)) {
    log_comp[, k] <- log(w[k]) + stats::dbeta(x, a[k], b[k], log = TRUE)
  }
  log_dens <- apply(log_comp, 1L, .geoess_log_sum_exp)
  if (log) log_dens else exp(log_dens)
}

#' @export
pmix <- function(mix, q, lower.tail = TRUE, log.p = FALSE) {
  UseMethod("pmix")
}

#' @method pmix betaMix
#' @export
pmix.betaMix <- function(mix, q, lower.tail = TRUE, log.p = FALSE) {
  .geoess_check_betaMix(mix)
  q <- as.numeric(q)
  w <- mix["w", ]
  a <- mix["a", ]
  b <- mix["b", ]
  K <- length(w)
  Nq <- length(q)

  log_comp <- matrix(NA_real_, nrow = Nq, ncol = K)
  for (k in seq_len(K)) {
    log_comp[, k] <- log(w[k]) + stats::pbeta(
      q, a[k], b[k], lower.tail = lower.tail, log.p = TRUE
    )
  }
  log_p <- apply(log_comp, 1L, .geoess_log_sum_exp)
  if (log.p) log_p else exp(log_p)
}

#' @export
qmix <- function(mix, p, lower.tail = TRUE, log.p = FALSE, tol = 1e-8,
                 maxiter = 100) {
  UseMethod("qmix")
}

#' @method qmix betaMix
#' @export
qmix.betaMix <- function(mix, p, lower.tail = TRUE, log.p = FALSE, tol = 1e-8,
                         maxiter = 100) {
  .geoess_check_betaMix(mix)
  p <- as.numeric(p)
  if (log.p) {
    if (any(p > 0)) stop("log.p is TRUE, so p must be <= 0")
  } else {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  res <- numeric(length(p))
  eps <- 1e-12
  lower <- eps
  upper <- 1 - eps

  for (i in seq_along(p)) {
    if (!log.p && p[i] <= 0) {
      res[i] <- if (lower.tail) 0 else 1
      next
    }
    if (!log.p && p[i] >= 1) {
      res[i] <- if (lower.tail) 1 else 0
      next
    }
    if (log.p && !is.finite(p[i]) && p[i] < 0) {
      res[i] <- if (lower.tail) 0 else 1
      next
    }
    if (log.p && p[i] == 0) {
      res[i] <- if (lower.tail) 1 else 0
      next
    }

    target <- p[i]
    f <- function(x) {
      pmix(mix, x, lower.tail = lower.tail, log.p = log.p) - target
    }

    f_lower <- f(lower)
    f_upper <- f(upper)
    if (!is.finite(f_lower)) f_lower <- -Inf
    if (!is.finite(f_upper)) f_upper <- Inf
    if (f_lower > 0 || f_upper < 0) {
      warning("Failed to bracket qmix; returning NA")
      res[i] <- NA_real_
      next
    }

    sol <- tryCatch(
      stats::uniroot(f, lower = lower, upper = upper, tol = tol,
                     maxiter = maxiter),
      error = function(e) NULL
    )
    if (is.null(sol)) {
      warning("qmix failed to converge; returning NA")
      res[i] <- NA_real_
    } else {
      res[i] <- sol$root
    }
  }

  res
}

#' @export
rmix <- function(mix, n) {
  UseMethod("rmix")
}

#' @method rmix betaMix
#' @export
rmix.betaMix <- function(mix, n) {
  .geoess_check_betaMix(mix)
  if (!is.numeric(n) || length(n) != 1L || n < 0) stop("n must be nonnegative")
  if (n == 0) return(numeric(0))

  w <- mix["w", ]
  a <- mix["a", ]
  b <- mix["b", ]
  K <- length(w)
  idx <- sample.int(K, size = n, replace = TRUE, prob = w)
  samp <- stats::rbeta(n, a[idx], b[idx])
  attr(samp, "component") <- idx
  samp
}

#' @method print betaMix
#' @export
print.betaMix <- function(x, ...) {
  cat("betaMix with", ncol(x), "components\n")
  print(unclass(x), ...)
  invisible(x)
}

.geoess_mixdist3 <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.matrix(args[[1]])) {
    mix <- args[[1]]
  } else {
    if (length(args) == 0L) stop("At least one component is required")
    lens <- vapply(args, length, integer(1))
    if (any(lens != 3L)) stop("All components must be length 3")
    mix <- do.call(cbind, args)
    if (is.null(colnames(mix))) {
      nm <- names(args)
      if (!is.null(nm) && all(nzchar(nm))) {
        colnames(mix) <- nm
      } else {
        colnames(mix) <- paste0("comp", seq_len(ncol(mix)))
      }
    }
  }

  if (!is.matrix(mix) || nrow(mix) != 3L) stop("mix must be a 3 x K matrix")
  if (any(!is.finite(mix))) stop("mix contains non-finite values")

  w <- mix[1, ]
  if (any(w < 0)) stop("weights must be nonnegative")
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("weights must sum to a positive value")
  mix[1, ] <- w / sw
  mix
}

.geoess_ms2beta <- function(m, s, drop = TRUE) {
  m <- as.numeric(m)
  s <- as.numeric(s)
  if (any(!is.finite(m)) || any(!is.finite(s))) stop("m and s must be finite")
  if (any(m <= 0 | m >= 1)) stop("m must be in (0, 1)")
  if (any(s <= 0)) stop("s must be positive")

  n <- m * (1 - m) / s^2 - 1
  if (any(n < 0)) stop("Invalid mean/sd: implies negative precision")

  ab <- cbind(a = n * m, b = n * (1 - m))
  if (drop) drop(ab) else ab
}

.geoess_mn2beta <- function(m, n, drop = TRUE) {
  m <- as.numeric(m)
  n <- as.numeric(n)
  if (any(!is.finite(m)) || any(!is.finite(n))) stop("m and n must be finite")
  if (any(m <= 0 | m >= 1)) stop("m must be in (0, 1)")
  if (any(n < 0)) stop("n must be nonnegative")

  ab <- cbind(a = n * m, b = n * (1 - m))
  if (drop) drop(ab) else ab
}

.geoess_check_betaMix <- function(mix) {
  if (!inherits(mix, "betaMix")) stop("mix must be a betaMix object")
  if (!is.matrix(mix) || nrow(mix) != 3L) stop("mix must be a 3 x K matrix")
  if (any(!is.finite(mix))) stop("mix contains non-finite values")
  invisible(TRUE)
}

#' Fit a Beta Mixture by EM
#'
#' Fits a mixture of Beta distributions to draws in (0, 1).
#'
#' @param draws Numeric vector or single-column matrix of draws.
#' @param K Integer. Number of mixture components.
#' @param init Character. Initialization method, "kmeans" or "random".
#' @param max_iter Integer. Maximum EM iterations.
#' @param tol Numeric. Convergence tolerance on log-likelihood.
#' @param min_weight Numeric. Minimum component weight.
#' @param min_shape Numeric. Minimum shape parameter to enforce positivity.
#' @param eps Numeric. Clamp draws to \code{[eps, 1-eps]} for stability.
#'
#' @return A \code{betaMix} object with fit diagnostics in attributes.
#' @export
geoess_mixfit_beta <- function(draws, K, init = c("kmeans", "random"),
                               max_iter = 500, tol = 1e-8,
                               min_weight = 1e-8, min_shape = 1e-4,
                               eps = 1e-8) {
  init <- match.arg(init)
  if (missing(K) || length(K) != 1L || K <= 0) stop("K must be a positive integer")

  x <- as.numeric(draws)
  if (!is.numeric(x)) stop("draws must be numeric")
  if (any(!is.finite(x))) stop("draws contain non-finite values")
  if (length(x) == 0L) stop("draws must be non-empty")

  x <- pmin(pmax(x, eps), 1 - eps)
  n <- length(x)
  if (n < K) stop("Number of draws must be >= K")

  if (!is.numeric(min_shape) || min_shape <= 0) stop("min_shape must be positive")
  if (!is.numeric(eps) || eps <= 0 || eps >= 0.5) stop("eps must be in (0, 0.5)")

  if (K == 1L) {
    m <- mean(x)
    v <- stats::var(x)
    ab <- .geoess_beta_ab_from_moments(m, v, min_shape = min_shape)
    a <- ab[1]
    b <- ab[2]
    logLik <- sum(stats::dbeta(x, a, b, log = TRUE))

    mix <- rbind(w = 1, a = a, b = b)
    class(mix) <- c("betaMix", "mix", "geoess_mixfit")
    attr(mix, "likelihood") <- "binomial"
    attr(mix, "family") <- "beta"
    attr(mix, "logLik") <- logLik
    attr(mix, "AIC") <- -2 * logLik + 2 * 2
    attr(mix, "BIC") <- -2 * logLik + log(n) * 2
    attr(mix, "iter") <- 1L
    attr(mix, "converged") <- TRUE
    attr(mix, "K") <- K
    return(mix)
  }

  if (init == "kmeans") {
    km <- stats::kmeans(x, centers = K, nstart = 5)
    cluster <- km$cluster
  } else {
    cluster <- sample(seq_len(K), n, replace = TRUE)
  }

  weights <- tabulate(cluster, nbins = K) / n
  a <- numeric(K)
  b <- numeric(K)
  for (k in seq_len(K)) {
    w_k <- as.numeric(cluster == k)
    if (sum(w_k) > 1) {
      m <- sum(w_k * x) / sum(w_k)
      v <- sum(w_k * (x - m)^2) / sum(w_k)
    } else {
      m <- mean(x)
      v <- stats::var(x)
    }
    ab <- .geoess_beta_ab_from_moments(m, v, min_shape = min_shape)
    a[k] <- ab[1]
    b[k] <- ab[2]
  }

  logLik <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    log_comp <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in seq_len(K)) {
      log_comp[, k] <- log(weights[k]) + stats::dbeta(x, a[k], b[k], log = TRUE)
    }
    log_mix <- apply(log_comp, 1L, .geoess_log_sum_exp)
    resp <- exp(log_comp - log_mix)
    logLik_new <- sum(log_mix)

    if (is.finite(logLik) &&
        abs(logLik_new - logLik) < tol * (1 + abs(logLik))) {
      converged <- TRUE
      logLik <- logLik_new
      break
    }
    logLik <- logLik_new

    Nk <- colSums(resp)
    weights <- pmax(Nk / n, min_weight)
    weights <- weights / sum(weights)

    for (k in seq_len(K)) {
      if (Nk[k] <= 1e-8) next

      w_k <- resp[, k]
      m <- sum(w_k * x) / Nk[k]
      v <- sum(w_k * (x - m)^2) / Nk[k]
      ab_init <- .geoess_beta_ab_from_moments(m, v, min_shape = min_shape)

      mle <- .geoess_beta_mle_weighted(
        x, w_k, ab_init[1], ab_init[2],
        tol = tol, max_iter = 100, min_shape = min_shape
      )
      a[k] <- mle$a
      b[k] <- mle$b
    }
  }

  mix <- rbind(w = weights, a = a, b = b)
  class(mix) <- c("betaMix", "mix", "geoess_mixfit")
  attr(mix, "likelihood") <- "binomial"
  attr(mix, "family") <- "beta"
  attr(mix, "logLik") <- logLik
  n_params <- (K - 1) + 2 * K
  attr(mix, "AIC") <- -2 * logLik + 2 * n_params
  attr(mix, "BIC") <- -2 * logLik + log(n) * n_params
  attr(mix, "iter") <- iter
  attr(mix, "converged") <- converged
  attr(mix, "K") <- K
  mix
}

#' Automated Beta Mixture Selection
#'
#' Fits beta mixtures across K and selects by AIC/BIC.
#'
#' @param draws Numeric vector or single-column matrix of draws.
#' @param K_grid Integer vector of candidate component counts.
#' @param criterion Character. "BIC" or "AIC".
#' @param ... Passed to \code{geoess_mixfit_beta}.
#'
#' @return A fitted \code{betaMix} object with selection metadata in attributes.
#' @export
geoess_automixfit_beta <- function(draws, K_grid = 1:6,
                                   criterion = c("BIC", "AIC"), ...) {
  criterion <- match.arg(criterion)
  if (!is.numeric(K_grid) || length(K_grid) == 0L) {
    stop("K_grid must be a non-empty numeric vector")
  }

  fits <- vector("list", length(K_grid))
  crit_vals <- rep(Inf, length(K_grid))

  for (i in seq_along(K_grid)) {
    K <- K_grid[i]
    fit <- tryCatch(
      geoess_mixfit_beta(draws, K = K, ...),
      error = function(e) NULL
    )
    fits[[i]] <- fit
    if (!is.null(fit)) {
      crit_vals[i] <- attr(fit, criterion)
    }
  }

  best_idx <- which.min(crit_vals)
  if (!is.finite(crit_vals[best_idx])) stop("All beta mixture fits failed")

  best_fit <- fits[[best_idx]]
  attr(best_fit, "criterion") <- criterion
  attr(best_fit, "K_selected") <- K_grid[best_idx]
  attr(best_fit, "criterion_values") <- crit_vals
  best_fit
}

.geoess_beta_ab_from_moments <- function(m, v, min_shape) {
  if (!is.finite(m) || !is.finite(v) || v <= 0) {
    return(c(min_shape, min_shape))
  }
  m <- min(max(m, 1e-6), 1 - 1e-6)
  phi <- m * (1 - m) / v - 1
  if (!is.finite(phi) || phi <= 0) phi <- min_shape
  a <- max(m * phi, min_shape)
  b <- max((1 - m) * phi, min_shape)
  c(a, b)
}

.geoess_beta_mle_weighted <- function(x, w, a_init, b_init, tol, max_iter,
                                      min_shape) {
  W <- sum(w)
  if (!is.finite(W) || W <= 0) {
    return(list(a = max(a_init, min_shape), b = max(b_init, min_shape),
                converged = FALSE))
  }

  logx <- log(x)
  log1x <- log(1 - x)
  S1 <- sum(w * logx)
  S2 <- sum(w * log1x)

  a <- max(a_init, min_shape)
  b <- max(b_init, min_shape)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    g1 <- S1 - W * (digamma(a) - digamma(a + b))
    g2 <- S2 - W * (digamma(b) - digamma(a + b))
    if (max(abs(g1), abs(g2)) < tol) {
      converged <- TRUE
      break
    }

    h11 <- -W * (trigamma(a) - trigamma(a + b))
    h22 <- -W * (trigamma(b) - trigamma(a + b))
    h12 <- W * trigamma(a + b)
    H <- matrix(c(h11, h12, h12, h22), nrow = 2)
    step <- tryCatch(solve(H, c(g1, g2)), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break

    step_scale <- 1
    repeat {
      a_new <- a - step_scale * step[1]
      b_new <- b - step_scale * step[2]
      if (is.finite(a_new) && is.finite(b_new) &&
          a_new > min_shape && b_new > min_shape) {
        break
      }
      step_scale <- step_scale / 2
      if (step_scale < 1e-6) {
        a_new <- max(a, min_shape)
        b_new <- max(b, min_shape)
        break
      }
    }
    a <- a_new
    b <- b_new
  }

  list(a = a, b = b, converged = converged)
}
