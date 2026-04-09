#' Fisher Information via Hessian
#'
#' Computes per-observation Fisher information from a log-likelihood Hessian.
#'
#' @param theta Numeric vector evaluation point.
#' @param loglik_fn Function returning scalar log-likelihood.
#' @param n_obs Integer number of observations.
#' @param symmetrize Logical. If TRUE, symmetrize the Hessian.
#' @param ... Passed to \code{loglik_fn}.
#'
#' @return Numeric matrix I1.
#' @export
geoess_I1_hessian <- function(theta, loglik_fn, n_obs, symmetrize = TRUE, ...) {
  if (!is.numeric(theta)) stop("theta must be numeric")
  if (!is.function(loglik_fn)) stop("loglik_fn must be a function")
  if (!is.numeric(n_obs) || length(n_obs) != 1L || n_obs <= 0) {
    stop("n_obs must be positive")
  }

  H <- numDeriv::hessian(func = loglik_fn, x = theta, ...)
  I1 <- -H / n_obs

  if (symmetrize) I1 <- .geoess_symmetrize(I1)
  if (!.geoess_check_pd(I1)) warning("I1 is not positive definite")
  I1
}

#' Fisher Information via Score Outer-Products
#'
#' Computes per-observation Fisher information from score vectors.
#'
#' @param theta Numeric vector evaluation point.
#' @param score_fn Function returning an n_obs x d score matrix.
#' @param n_obs Integer number of observations.
#' @param ... Passed to \code{score_fn}.
#'
#' @return Numeric matrix I1.
#' @export
geoess_I1_score <- function(theta, score_fn, n_obs, ...) {
  if (!is.numeric(theta)) stop("theta must be numeric")
  if (!is.function(score_fn)) stop("score_fn must be a function")
  if (!is.numeric(n_obs) || length(n_obs) != 1L || n_obs <= 0) {
    stop("n_obs must be positive")
  }

  scores <- score_fn(theta, ...)
  scores <- as.matrix(scores)
  if (!is.numeric(scores)) stop("scores must be numeric")
  if (nrow(scores) != n_obs) stop("scores must have n_obs rows")

  I1 <- crossprod(scores) / n_obs
  I1 <- .geoess_symmetrize(I1)
  if (!.geoess_check_pd(I1)) warning("I1 is not positive definite")
  I1
}

geoess_I1_from_fn <- function(theta, I1_fn, symmetrize = TRUE) {
  if (!is.function(I1_fn)) stop("I1_fn must be a function")
  I1 <- I1_fn(theta)
  I1 <- .geoess_as_matrix(I1, "I1")
  .geoess_check_square(I1, "I1")
  if (symmetrize) I1 <- .geoess_symmetrize(I1)
  I1
}

geoess_I1_from_hessian <- function(theta, loglik_fn, n_obs, ...) {
  geoess_I1_hessian(theta, loglik_fn, n_obs, ...)
}
