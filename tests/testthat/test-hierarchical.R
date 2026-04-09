ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "hier.R"))
}

test_that("schur_theta_given_phi matches manual Schur complement", {
  Htt <- matrix(2, 1, 1)
  Htp <- matrix(1, 1, 1)
  Hpp <- matrix(3, 1, 1)
  res <- schur_theta_given_phi(Htt, Htp, Hpp)
  expect_equal(res$Ipi_eff, matrix(5 / 3, 1, 1))
  expect_equal(res$Ipi_between, Htt)
  expect_equal(res$Ipi_borrowed, matrix(-1 / 3, 1, 1))
})

test_that("schur complement returns symmetric blocks", {
  set.seed(1)
  Htt <- crossprod(matrix(rnorm(9), 3, 3))
  Htp <- matrix(rnorm(6), 3, 2)
  Hpp <- crossprod(matrix(rnorm(4), 2, 2)) + diag(0.5, 2)
  res <- schur_theta_given_phi(Htt, Htp, Hpp)
  expect_equal(res$Ipi_eff, t(res$Ipi_eff))
  expect_equal(res$Ipi_between, t(res$Ipi_between))
  expect_equal(res$Ipi_borrowed, t(res$Ipi_borrowed))
})

test_that("Normal-Normal hierarchical ESS decomposition holds", {
  n <- 5
  sigma <- 2
  tau <- 1.5
  s0 <- 3

  Htt <- diag(1 / tau^2, n)
  Htp <- matrix(-1 / tau^2, n, 1)
  Hpp <- matrix(n / tau^2 + 1 / s0^2, 1, 1)

  I1 <- diag(1 / sigma^2, n)
  blocks <- geoess_Ipi_hier(Htt, Htp, Hpp, return_components = TRUE)

  res <- geoess_ess_hier(I1 = I1, Ipi_blocks = blocks,
                         theta_ref = rep(0, n), phi_ref = 0,
                         method = "trace")

  expect_true(res$ess_borrowed < 0)
  expect_equal(res$ess_total, res$ess_between + res$ess_borrowed, tolerance = 1e-10)
})

test_that("Directional ESS shows anisotropy in hierarchical model", {
  n <- 4
  sigma <- 2
  tau <- 1
  s0 <- 2

  Htt <- diag(1 / tau^2, n)
  Htp <- matrix(-1 / tau^2, n, 1)
  Hpp <- matrix(n / tau^2 + 1 / s0^2, 1, 1)

  I1 <- diag(1 / sigma^2, n)
  blocks <- geoess_Ipi_hier(Htt, Htp, Hpp, return_components = TRUE)

  v1 <- rep(1, n)
  v1 <- v1 / sqrt(sum(v1^2))
  v2 <- c(1, rep(0, n - 1))

  res1 <- geoess_ess_hier(I1 = I1, Ipi_blocks = blocks,
                          theta_ref = rep(0, n), phi_ref = 0,
                          method = "dir", v = v1)
  res2 <- geoess_ess_hier(I1 = I1, Ipi_blocks = blocks,
                          theta_ref = rep(0, n), phi_ref = 0,
                          method = "dir", v = v2)

  expect_true(abs(res1$ess_total - res2$ess_total) > 1e-8)
})
