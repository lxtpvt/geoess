#' Predictive Effective Sample Size (Template)
#'
#' Matches a target KL divergence on a grid of sample sizes.
#'
#' @param target_kl Numeric. Target KL divergence between predictive distributions.
#' @param n_grid Numeric vector. Candidate sample sizes.
#' @param kl_calc_fn Function. Function(n) returning KL divergence for given n.
#'
#' @return A list with the best-matching n, grid values, and KL diagnostics.
#'
#' @examples
#' kl_fn <- function(n) 1 / (n + 1)
#' predictive_ess_template(0.2, n_grid = 1:10, kl_calc_fn = kl_fn)$n_eff
#' @export
predictive_ess_template <- function(target_kl, n_grid, kl_calc_fn) {
  if (!is.numeric(target_kl) || length(target_kl) != 1L) {
    stop("target_kl must be a numeric scalar")
  }
  if (!is.numeric(n_grid) || length(n_grid) == 0L) {
    stop("n_grid must be a numeric vector")
  }
  if (!is.function(kl_calc_fn)) stop("kl_calc_fn must be a function")

  kl_values <- vapply(n_grid, kl_calc_fn, numeric(1))
  diffs <- abs(kl_values - target_kl)
  best_idx <- which.min(diffs)

  list(
    n_eff = n_grid[best_idx],
    grid = n_grid,
    kl_values = kl_values,
    min_diff = diffs[best_idx]
  )
}

#' Predictive ESS for Normal Mean with Known Variance
#'
#' Matches predictive KL divergences for a Normal mean model using analytic KL.
#'
#' @param sigma Numeric. Data standard deviation.
#' @param s0 Numeric. Prior standard deviation.
#' @param n0 Numeric. Baseline sample size for predictive comparison.
#' @param tol Numeric. Root-finding tolerance.
#' @param max_n Numeric. Maximum n for bracketing.
#'
#' @return A list containing n_eff, target KL, matched KL, and max_kl.
#'
#' @details
#' The KL calculation assumes equal predictive means (e.g., evaluated at the prior mean),
#' so the divergence depends only on predictive variances. Matching can be infeasible
#' when the target KL exceeds the maximum achievable KL for the baseline \code{n0};
#' in that case the function returns \code{NA} with a warning.
#'
#' @examples
#' predictive_ess_normal_mu(sigma = 2, s0 = 2, n0 = 1)
#' @export
predictive_ess_normal_mu <- function(sigma, s0, n0 = 1, tol = 1e-8, max_n = 1e6) {
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("sigma must be positive")
  }
  if (!is.numeric(s0) || length(s0) != 1L || s0 <= 0) {
    stop("s0 must be positive")
  }
  if (!is.numeric(n0) || length(n0) != 1L || n0 <= 0) {
    stop("n0 must be positive")
  }

  s_n2 <- function(n) 1 / (1 / s0^2 + n / sigma^2)
  v_pred <- function(n) sigma^2 + s_n2(n)

  v0 <- sigma^2 + s0^2
  v_n0 <- v_pred(n0)
  target_kl <- .geoess_kl_normal_var(v0, v_n0)
  max_kl <- .geoess_kl_normal_var(v_n0, sigma^2)

  if (target_kl > max_kl + tol) {
    warning("Target KL exceeds maximum achievable KL for this n0; returning NA")
    return(list(
      n_eff = NA_real_,
      target_kl = target_kl,
      matched_kl = NA_real_,
      max_kl = max_kl,
      n0 = n0,
      sigma = sigma,
      s0 = s0
    ))
  }

  f <- function(n_eff) {
    .geoess_kl_normal_var(v_n0, v_pred(n0 + n_eff)) - target_kl
  }

  upper <- max(1, n0)
  f_upper <- f(upper)
  while (f_upper < 0 && upper < max_n) {
    upper <- upper * 2
    f_upper <- f(upper)
  }

  if (f_upper < 0) {
    warning("Unable to bracket predictive ESS; returning NA")
    return(list(
      n_eff = NA_real_,
      target_kl = target_kl,
      matched_kl = NA_real_,
      max_kl = max_kl,
      n0 = n0,
      sigma = sigma,
      s0 = s0
    ))
  }

  sol <- uniroot(f, lower = 0, upper = upper, tol = tol)
  n_eff <- sol$root
  matched_kl <- .geoess_kl_normal_var(v_n0, v_pred(n0 + n_eff))

  list(
    n_eff = n_eff,
    target_kl = target_kl,
    matched_kl = matched_kl,
    max_kl = max_kl,
    n0 = n0,
    sigma = sigma,
    s0 = s0
  )
}

#' Predictive ESS on a Discrete Grid
#'
#' Matches KL divergence between a baseline and target predictive distribution
#' to the closest candidate distribution on a grid.
#'
#' @param target_p Numeric vector. Target predictive probabilities.
#' @param candidate_grid List or matrix of candidate predictive distributions.
#' @param baseline_p Numeric vector. Baseline predictive probabilities.
#' @param n_grid Optional numeric vector indexing the candidates.
#' @param eps Numeric. Small value to stabilize zero probabilities.
#'
#' @return A list with the matched n_eff, target KL, and grid KL values.
#'
#' @examples
#' baseline <- c(0.5, 0.5)
#' target <- c(0.7, 0.3)
#' candidates <- rbind(c(0.55, 0.45), c(0.6, 0.4), c(0.7, 0.3))
#' predictive_ess_discrete(target, candidates, baseline, n_grid = 1:3)$n_eff
#' @export
predictive_ess_discrete <- function(target_p, candidate_grid, baseline_p,
                                    n_grid = NULL, eps = 1e-12) {
  if (missing(baseline_p)) stop("baseline_p is required")

  target_p <- as.numeric(target_p)
  baseline_p <- as.numeric(baseline_p)

  if (length(target_p) != length(baseline_p)) {
    stop("target_p and baseline_p must be the same length")
  }

  candidates <- NULL
  if (is.list(candidate_grid)) {
    candidates <- candidate_grid
  } else if (is.matrix(candidate_grid)) {
    if (ncol(candidate_grid) == length(target_p)) {
      candidates <- lapply(seq_len(nrow(candidate_grid)), function(i) candidate_grid[i, ])
    } else if (nrow(candidate_grid) == length(target_p)) {
      candidates <- lapply(seq_len(ncol(candidate_grid)), function(i) candidate_grid[, i])
    } else {
      stop("candidate_grid has incompatible dimensions")
    }
  } else {
    stop("candidate_grid must be a list or matrix")
  }

  n_candidates <- length(candidates)
  if (is.null(n_grid)) n_grid <- seq_len(n_candidates)
  if (length(n_grid) != n_candidates) stop("n_grid length mismatch")

  target_kl <- .geoess_kl_discrete(baseline_p, target_p, eps = eps)
  grid_kl <- vapply(candidates, function(p) {
    .geoess_kl_discrete(baseline_p, p, eps = eps)
  }, numeric(1))

  diffs <- abs(grid_kl - target_kl)
  best_idx <- which.min(diffs)

  list(
    n_eff = n_grid[best_idx],
    target_kl = target_kl,
    grid_kl = grid_kl,
    best_idx = best_idx,
    n_grid = n_grid
  )
}
