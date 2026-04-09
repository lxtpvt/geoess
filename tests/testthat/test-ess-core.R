ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {

  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "models_normal.R"))
}

test_that("geoess_normal_mu matches analytic value", {
  expect_equal(geoess_normal_mu(2, 4), 0.25)
})

test_that("geoess_directional matches ratio definition", {
  I1 <- diag(c(2, 3))
  Ipi <- diag(c(1, 6))
  v <- c(2, -1)
  expected <- as.numeric(t(v) %*% Ipi %*% v / (t(v) %*% I1 %*% v))
  expect_equal(geoess_directional(I1, Ipi, v), expected)
})

test_that("geoess_eigen nonnegative for PSD Ipi and PD I1", {
  I1 <- diag(c(2, 4))
  Ipi <- diag(2)
  eig <- geoess_eigen(I1, Ipi)
  expect_true(all(eig > -1e-8))
})

test_that("geoess_curvature matches manual trace", {
  I1 <- diag(2)
  Ipi <- diag(c(0.5, 0.5))
  res <- geoess_curvature(I1, Ipi)
  expect_equal(res$ess, 0.5)
})

test_that("geoess_trace matches geoess_curvature", {
  I1 <- diag(c(2, 3))
  Ipi <- diag(c(1, 1))
  expect_equal(geoess_trace(I1, Ipi), geoess_curvature(I1, Ipi)$ess)
})

test_that("geoess_ess returns weighted component ESS for mixtures", {
  I1 <- diag(c(2, 2))
  mix <- structure(
    list(
      weights = c(0.25, 0.75),
      mu = list(c(0, 0), c(0, 0)),
      Sigma = list(diag(2), diag(c(2, 2)))
    ),
    class = "geoess_mix"
  )
  res <- geoess_ess(I1 = I1, Ipi = mix, method = "trace")
  expect_equal(res$ess_components, c(0.5, 0.25), tolerance = 1e-12)
  expect_equal(res$ess_weighted, 0.3125, tolerance = 1e-12)
})
