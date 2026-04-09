test_that("Schur complement matches manual formula", {
  set.seed(1)
  A <- matrix(rnorm(16), 4, 4)
  M <- crossprod(A) + diag(0.1, 4)
  idx_a <- c(1, 2)
  idx_b <- c(3, 4)

  sch <- geoess:::.geoess_schur(M, idx_a, idx_b)
  manual <- M[idx_a, idx_a] - M[idx_a, idx_b] %*% solve(M[idx_b, idx_b], M[idx_b, idx_a])

  expect_equal(sch, manual, tolerance = 1e-8)
  expect_equal(sch, t(sch), tolerance = 1e-8)
})

test_that("Conditional ESS reduces to marginal ratio for diagonal matrices", {
  I1 <- diag(c(2, 3))
  Ipi <- diag(c(4, 9))
  res <- geoess_conditional(I1, Ipi, idx_a = 1, idx_b = 2)
  expect_equal(as.numeric(res$ess_cond), 4 / 2, tolerance = 1e-8)
})

test_that("Conditional ESS differs from marginal under correlation", {
  I1 <- diag(c(1, 1))
  Ipi <- matrix(c(2, 1, 1, 2), 2, 2)
  res <- geoess_conditional(I1, Ipi, idx_a = 1, idx_b = 2)
  expect_true(is.finite(res$ess_cond))
  expect_true(abs(res$ess_cond - (Ipi[1, 1] / I1[1, 1])) > 1e-6)
})

test_that("Spectrum matches diagonal ratios and reconstructs M", {
  I1 <- diag(c(1, 2))
  Ipi <- diag(c(4, 6))
  spec <- geoess_spectrum(I1, Ipi)
  expect_equal(sort(spec$values), sort(c(4, 3)), tolerance = 1e-8)
  recon <- spec$vectors %*% diag(spec$values) %*% t(spec$vectors)
  expect_equal(recon, spec$M, tolerance = 1e-8)
})

test_that("Global conditional ESS averages local values", {
  set.seed(1)
  draws <- matrix(rnorm(200), ncol = 2)
  I1 <- diag(c(2, 5))
  Ipi <- diag(c(4, 10))
  res <- geoess_conditional_global(draws, I1, Ipi, idx_a = 1, idx_b = 2, nsim = 200)
  expect_equal(as.numeric(res$estimate), 4 / 2, tolerance = 1e-8)
})
