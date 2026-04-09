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

test_that("demo_normal_mu returns matching ESS", {
  res <- demo_normal_mu(n = 50, mu_true = 1, sigma = 2, m0 = 0, s0 = 4, seed = 123)
  expect_true(res$match)
  expect_equal(res$ess_geometric, 0.25, tolerance = 1e-12)
  expect_equal(res$ess_theoretical, 0.25, tolerance = 1e-12)
})
