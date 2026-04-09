.geoess_log_sum_exp <- function(x) {
  x <- as.numeric(x)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.geoess_logsumexp <- function(x) {
  .geoess_log_sum_exp(x)
}

.geoess_softmax_log <- function(log_w) {
  log_w <- as.numeric(log_w)
  lse <- .geoess_log_sum_exp(log_w)
  exp(log_w - lse)
}

.geoess_check_pd <- function(A, tol = 1e-8) {
  if (!is.matrix(A)) A <- as.matrix(A)
  if (nrow(A) == 0 || ncol(A) == 0) return(FALSE)
  if (nrow(A) != ncol(A)) return(FALSE)
  ev <- eigen(.geoess_symmetrize(A), symmetric = TRUE, only.values = TRUE)$values
  all(ev > -tol)
}

.geoess_log_mvn <- function(x, mu, Sigma, chol = NULL) {
  x <- as.numeric(x)
  mu <- as.numeric(mu)
  if (length(x) != length(mu)) stop("x and mu must have same length")

  if (!is.null(chol)) {
    R <- chol
    if (!is.matrix(R) || nrow(R) != ncol(R) || nrow(R) != length(x)) {
      stop("chol has incompatible dimensions")
    }
    z <- backsolve(R, x - mu, transpose = TRUE)
    quad <- sum(z^2)
    logdet <- 2 * sum(log(diag(R)))
    return(-0.5 * (length(x) * log(2 * pi) + logdet + quad))
  }

  mvtnorm::dmvnorm(x, mean = mu, sigma = Sigma, log = TRUE)
}
