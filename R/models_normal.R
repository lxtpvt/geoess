#' Fisher Information for Normal Mean (Per Observation)
#'
#' @param sigma Numeric. Data standard deviation.
#'
#' @return 1x1 matrix for \eqn{1/\sigma^2}.
#'
#' @examples
#' I1_normal_mu(2)
#' @export
I1_normal_mu <- function(sigma) {
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("sigma must be a positive scalar")
  }
  matrix(1 / sigma^2, nrow = 1L, ncol = 1L)
}

#' Prior Curvature for Normal Mean (Normal Prior)
#'
#' @param s0 Numeric. Prior standard deviation.
#'
#' @return 1x1 matrix for \eqn{1/s_0^2}.
#'
#' @examples
#' Ipi_normal_mu(4)
#' @export
Ipi_normal_mu <- function(s0) {
  if (!is.numeric(s0) || length(s0) != 1L || s0 <= 0) {
    stop("s0 must be a positive scalar")
  }
  matrix(1 / s0^2, nrow = 1L, ncol = 1L)
}

#' Geometric ESS for Normal Mean with Known Variance
#'
#' @param sigma Numeric. Data standard deviation.
#' @param s0 Numeric. Prior standard deviation.
#'
#' @return Numeric scalar ESS, \eqn{\sigma^2/s_0^2}.
#'
#' @examples
#' geoess_normal_mu(sigma = 2, s0 = 4)
#' @export
geoess_normal_mu <- function(sigma, s0) {
  I1 <- I1_normal_mu(sigma)
  Ipi <- Ipi_normal_mu(s0)
  geoess_curvature(I1, Ipi)$ess
}

#' Demonstrate Geometric ESS for Normal Mean
#'
#' Simulates data and compares geometric ESS with the analytic value.
#'
#' @param n Integer. Sample size.
#' @param mu_true Numeric. True mean used to simulate data.
#' @param sigma Numeric. Data standard deviation.
#' @param m0 Numeric. Prior mean.
#' @param s0 Numeric. Prior standard deviation.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list with ESS values, posterior summaries, and a match flag.
#'
#' @examples
#' demo_normal_mu(n = 50, mu_true = 1, sigma = 2, m0 = 0, s0 = 4, seed = 123)
#' @export
demo_normal_mu <- function(n = 100, mu_true = 0, sigma = 1, m0 = 0, s0 = 1,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(n) || length(n) != 1L || n <= 0) stop("n must be positive")
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("sigma must be positive")
  }
  if (!is.numeric(s0) || length(s0) != 1L || s0 <= 0) {
    stop("s0 must be positive")
  }

  y <- rnorm(n, mean = mu_true, sd = sigma)
  xbar <- mean(y)

  I1 <- I1_normal_mu(sigma)
  Ipi <- Ipi_normal_mu(s0)
  ess_geometric <- geoess_curvature(I1, Ipi)$ess
  ess_theoretical <- sigma^2 / s0^2
  match <- abs(ess_geometric - ess_theoretical) < 1e-12

  post_prec <- n / sigma^2 + 1 / s0^2
  post_mean <- (n * xbar / sigma^2 + m0 / s0^2) / post_prec
  post_var <- 1 / post_prec

  list(
    ess_geometric = ess_geometric,
    ess_theoretical = ess_theoretical,
    match = match,
    xbar = xbar,
    post_mean = post_mean,
    post_var = post_var
  )
}
