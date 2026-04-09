ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "predictive.R"))
}

test_that("predictive_ess_normal_mu matches KL target", {
  res <- predictive_ess_normal_mu(sigma = 2, s0 = 2, n0 = 1)
  expect_true(is.finite(res$n_eff))
  expect_true(res$n_eff > 0)
  expect_true(abs(res$matched_kl - res$target_kl) < 1e-6)
})

test_that("predictive_ess_discrete selects closest candidate", {
  baseline <- c(0.5, 0.5)
  target <- c(0.7, 0.3)
  candidates <- rbind(c(0.55, 0.45), c(0.6, 0.4), c(0.7, 0.3))
  res <- predictive_ess_discrete(target, candidates, baseline, n_grid = 1:3)
  expect_equal(res$n_eff, 3)
})
