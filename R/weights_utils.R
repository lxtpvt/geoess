#' Weight Utilities for Global ESS
#'
#' Internal helpers for normalizing weights and computing diagnostics.
#'
#' @keywords internal
.geoess_normalize_weights <- function(w, n) {
  if (is.null(w)) return(rep(1 / n, n))
  w <- as.numeric(w)
  if (length(w) != n) stop("weights length does not match number of draws")
  if (any(!is.finite(w))) stop("weights contain non-finite values")
  if (any(w < 0)) stop("weights must be nonnegative")
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("weights must sum to a positive value")
  w / sw
}

.geoess_normalize_logw <- function(logw) {
  logw <- as.numeric(logw)
  lse <- .geoess_log_sum_exp(logw)
  w <- exp(logw - lse)
  list(weights = w, log_norm = lse)
}

.geoess_is_ess <- function(w) {
  w <- as.numeric(w)
  w <- w / sum(w)
  1 / sum(w^2)
}

.geoess_weight_summary <- function(w) {
  w <- as.numeric(w)
  w <- w / sum(w)
  list(
    ess = .geoess_is_ess(w),
    min = min(w),
    max = max(w),
    cv = stats::sd(w) / mean(w)
  )
}

.geoess_weighted_mean <- function(x, w) {
  sum(w * x)
}

.geoess_weighted_var <- function(x, w) {
  mu <- .geoess_weighted_mean(x, w)
  sum(w * (x - mu)^2)
}

.geoess_weighted_quantile <- function(x, w, probs = c(0.05, 0.5, 0.95)) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  if (length(x) != length(w)) stop("x and w must have same length")
  if (any(w < 0)) stop("weights must be nonnegative")
  w <- w / sum(w)

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w)

  vapply(probs, function(p) {
    idx <- which(cw >= p)[1]
    x[idx]
  }, numeric(1))
}
