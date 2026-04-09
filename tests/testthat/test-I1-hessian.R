ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "fisher_info.R"))
}

test_that("geoess_I1_hessian matches Normal mean Fisher", {
  set.seed(1)
  y <- rnorm(10, mean = 0, sd = 2)
  loglik_fn <- function(theta, y, sigma) {
    sum(dnorm(y, mean = theta, sd = sigma, log = TRUE))
  }
  I1 <- geoess_I1_hessian(0, loglik_fn, n_obs = length(y), y = y, sigma = 2)
  expect_equal(as.numeric(I1), 1 / 4, tolerance = 1e-4)
})
