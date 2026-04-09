test_that("geoess_hier_global matches Normal-Normal global ESS constants", {
  set.seed(123)
  K <- 4
  sigma <- 2
  tau <- 0.5
  s0 <- 1.0
  n_draws <- 300

  draws <- cbind(matrix(rnorm(n_draws * K), ncol = K), rnorm(n_draws))
  pars_theta <- paste0("theta[", seq_len(K), "]")
  pars_phi <- "mu"
  colnames(draws) <- c(pars_theta, pars_phi)

  Ipi_fn <- function(Theta, ...) {
    Htt <- diag(1 / tau^2, K)
    Htp <- matrix(-1 / tau^2, nrow = K, ncol = 1)
    Hpp <- matrix(K / tau^2 + 1 / s0^2, nrow = 1, ncol = 1)
    list(
      Ipi_thetatheta = Htt,
      Ipi_thetaphi = Htp,
      Ipi_phiphi = Hpp
    )
  }
  I1_theta_fn <- function(theta_i, phi, ...) matrix(1 / sigma^2, 1, 1)
  I1_phi_fn <- function(phi, theta, ...) matrix(1 / tau^2, 1, 1)

  res <- geoess_hier_global(
    Q = "custom",
    draws = draws,
    pars_theta = pars_theta,
    pars_phi = pars_phi,
    group_index = seq_len(K),
    Ipi_provider = "analytic",
    Ipi_fn = Ipi_fn,
    I1_provider = "analytic",
    I1_theta_fn = I1_theta_fn,
    I1_phi_fn = I1_phi_fn
  )

  target_phi <- tau^2 / s0^2 + K
  target_theta <- (sigma^2 / tau^2) * (1 - 1 / (K + tau^2 / s0^2))

  expect_equal(res$ess_phi_global$estimate, target_phi, tolerance = 1e-8)
  expect_equal(as.numeric(res$ess_theta_cond_global),
               rep(target_theta, K),
               tolerance = 1e-6)

  w <- runif(n_draws)
  res_w <- geoess_hier_global(
    Q = "custom",
    draws = draws,
    weights = w,
    pars_theta = pars_theta,
    pars_phi = pars_phi,
    group_index = seq_len(K),
    Ipi_provider = "analytic",
    Ipi_fn = Ipi_fn,
    I1_provider = "analytic",
    I1_theta_fn = I1_theta_fn,
    I1_phi_fn = I1_phi_fn
  )

  expect_equal(res_w$ess_phi_global$estimate, target_phi, tolerance = 1e-8)
  expect_equal(as.numeric(res_w$ess_theta_cond_global),
               rep(target_theta, K),
               tolerance = 1e-6)
})
