ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "beta_binomial_canonical.R"))
}

test_that("geoess_beta_binomial_canonical returns classical ESS a+b", {
  expect_equal(geoess_beta_binomial_canonical(a = 3, b = 7), 10)
  expect_equal(geoess_beta_binomial_canonical(a = 1, b = 1), 2)
})

test_that("ratio_eta is constant on the canonical scale", {
  det <- geoess_beta_binomial_canonical(a = 3, b = 7, return = "details")
  grid <- c(-5, -1, 0, 1, 5)

  expect_equal(det$ratio_eta(grid), rep(10, length(grid)), tolerance = 1e-12)
})

test_that("Ipi_eta divided by I1_eta equals a+b numerically", {
  a <- 4
  b <- 9
  det <- geoess_beta_binomial_canonical(a = a, b = b, return = "details")
  grid <- c(-5, -1, 0, 1, 5)

  ratio <- det$Ipi_eta(grid) / det$I1_eta(grid)
  expect_equal(ratio, rep(a + b, length(grid)), tolerance = 1e-12)
})

test_that("ESS remains per-trial when n is supplied", {
  det <- geoess_beta_binomial_canonical(a = 2, b = 5, n = 40, return = "details")
  grid <- c(-2, 0, 2)

  expect_equal(det$ess_local, 7, tolerance = 1e-12)
  expect_equal(det$ess_global, 7, tolerance = 1e-12)
  expect_equal(det$ratio_eta(grid), rep(7, length(grid)), tolerance = 1e-12)
  expect_true(grepl("Per-trial ESS", det$notes, fixed = TRUE))
  expect_true(grepl("n = 40", det$notes, fixed = TRUE))
})
