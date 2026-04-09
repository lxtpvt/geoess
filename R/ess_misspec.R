#' Misspecified Effective Sample Size using Godambe Information
#'
#' Computes ESS under model misspecification using the Godambe metric.
#'
#' @param I Numeric matrix. Fisher curvature (expected negative Hessian).
#' @param J Numeric matrix. Variance of the score function.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param symmetrize Logical. If TRUE, symmetrize inputs before calculations.
#'
#' @return A list with scalar, eigen, and directional ESS as in \code{geoess_curvature}.
#'
#' @examples
#' I <- diag(c(2, 3))
#' J <- I
#' Ipi <- diag(c(1, 1))
#' geoess_godambe(I, J, Ipi)$ess
#' @export
geoess_godambe <- function(I, J, Ipi, symmetrize = TRUE) {
  I <- .geoess_as_matrix(I, "I")
  J <- .geoess_as_matrix(J, "J")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I, "I")
  .geoess_check_square(J, "J")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I, J, "I", "J")
  .geoess_check_same_dim(I, Ipi, "I", "Ipi")

  if (symmetrize) {
    I <- .geoess_symmetrize(I)
    J <- .geoess_symmetrize(J)
    Ipi <- .geoess_symmetrize(Ipi)
  }

  tmp <- .geoess_solve(J, I)
  G <- I %*% tmp
  geoess_curvature(G, Ipi, symmetrize = symmetrize)
}

.geoess_I1_from_loglik <- function(loglik_fn, theta, n = 1, ...) {
  if (!is.function(loglik_fn)) stop("loglik_fn must be a function")
  if (!is.numeric(theta)) stop("theta must be numeric")
  if (length(n) != 1L || n <= 0) stop("n must be positive")

  H <- numDeriv::hessian(func = loglik_fn, x = theta, ...)
  I1 <- -H / n
  .geoess_symmetrize(I1)
}

.geoess_J_from_scores <- function(scores) {
  S <- as.matrix(scores)
  if (!is.numeric(S)) stop("scores must be numeric")
  if (nrow(S) < 2L) stop("scores must have at least 2 rows")

  centered <- sweep(S, 2, colMeans(S), FUN = "-")
  crossprod(centered) / (nrow(S) - 1)
}
