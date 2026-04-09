ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "mix_beta.R"))
}

test_that("geoess_mixfit_beta fits a single beta component", {
  set.seed(123)
  draws <- rbeta(300, 2, 6)
  fit <- geoess_mixfit_beta(draws, K = 1)
  expect_true(inherits(fit, "betaMix"))
  expect_equal(sum(fit["w", ]), 1, tolerance = 1e-12)
  expect_true(all(fit["a", ] > 0))
  expect_true(all(fit["b", ] > 0))
  expect_true(is.finite(attr(fit, "logLik")))
})

test_that("geoess_automixfit_beta returns selection metadata", {
  set.seed(1)
  draws <- rbeta(200, 3, 5)
  fit <- geoess_automixfit_beta(draws, K_grid = 1:2, criterion = "BIC")
  expect_true(inherits(fit, "betaMix"))
  expect_true(is.numeric(attr(fit, "K_selected")))
  expect_true(length(attr(fit, "criterion_values")) == 2)
})
