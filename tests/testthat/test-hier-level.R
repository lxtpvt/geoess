test_that("Normal–Normal hierarchical level ESS matches analytic formulas", {
  K <- 5
  sigma <- 2
  tau <- 0.5
  s0 <- 1.2

  Iunit_theta <- matrix(1 / sigma^2, 1, 1)
  Ipi_theta <- matrix(1 / tau^2, 1, 1)
  ess_theta <- geoess_block_trace(Iunit_theta, Ipi_theta)
  expect_equal(as.numeric(ess_theta), sigma^2 / tau^2, tolerance = 1e-8)

  Iunit_phi <- matrix(1 / tau^2, 1, 1)
  Ipi_phi <- matrix(K / tau^2 + 1 / s0^2, 1, 1)
  ess_phi <- geoess_block_trace(Iunit_phi, Ipi_phi)
  expect_equal(as.numeric(ess_phi), tau^2 / s0^2 + K, tolerance = 1e-8)
})

test_that("Conditional ESS via Schur complement matches analytic formula", {
  K <- 6
  sigma <- 1.5
  tau <- 0.4
  s0 <- 0.8

  Iunit_full <- diag(c(rep(1 / sigma^2, K), 1 / tau^2))

  Htt <- diag(1 / tau^2, K)
  Htp <- matrix(-1 / tau^2, nrow = K, ncol = 1)
  Hpp <- matrix(K / tau^2 + 1 / s0^2, nrow = 1, ncol = 1)
  Ipi_full <- rbind(
    cbind(Htt, Htp),
    cbind(t(Htp), Hpp)
  )

  idx_theta1 <- 1
  idx_mu <- K + 1
  cond <- geoess_conditional(Iunit_full, Ipi_full, idx_a = idx_theta1, idx_b = idx_mu)

  target <- (sigma^2 / tau^2) * (1 - 1 / (K + tau^2 / s0^2))
  expect_equal(as.numeric(cond$ess_cond), target, tolerance = 1e-8)
})

test_that("Spectrum for diagonal blocks matches coordinate ratios", {
  Iunit <- diag(c(2, 5))
  Ipi <- diag(c(6, 10))
  spec <- geoess_block_spectrum(Iunit, Ipi)
  expect_equal(sort(spec$values), sort(c(3, 2)), tolerance = 1e-8)
})
