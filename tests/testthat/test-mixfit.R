ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "mixfit_mvnorm_em.R"))
}

test_that("geoess_mixfit returns a valid mixture", {
  set.seed(123)
  draws <- cbind(c(rnorm(100, mean = -1), rnorm(100, mean = 2)))
  fit <- geoess_mixfit(draws, family = "norm", K = 2, init = "kmeans", max_iter = 200)
  expect_true(inherits(fit, "geoess_mix"))
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  expect_equal(length(fit$mu), 2)
  expect_true(is.finite(fit$logLik))
})

test_that("geoess_automixfit prefers K=1 for unimodal draws", {
  set.seed(42)
  draws <- cbind(rnorm(500))
  fit <- geoess_automixfit(draws, family = "norm", K_grid = 1:2, criterion = "BIC")
  expect_equal(fit$K_selected, 1)
})
