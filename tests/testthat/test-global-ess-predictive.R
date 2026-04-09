ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "weights_utils.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "models_normal.R"))
  source(file.path(pkg_root, "R", "global_ess.R"))
  source(file.path(pkg_root, "R", "global_ess_result.R"))
}

test_that("geoess_global_predictive returns finite ESS with weights", {
  sigma <- 1.5
  s0 <- 2.5
  I1 <- I1_normal_mu(sigma)
  Ipi <- Ipi_normal_mu(s0)

  set.seed(10)
  draws <- matrix(rnorm(1000, mean = 0, sd = s0), ncol = 1)
  yrep <- 0.2

  loglik_point <- function(theta, yrep) {
    theta <- as.numeric(theta)
    dnorm(yrep, mean = theta, sd = sigma, log = TRUE)
  }

  res <- geoess_global_predictive(
    prior = draws,
    I1 = I1,
    yrep = yrep,
    loglik_point = loglik_point,
    ipi_source = "fixed",
    Ipi = Ipi
  )

  expect_true(is.finite(res$estimate))
  expect_equal(res$estimate, sigma^2 / s0^2, tolerance = 1e-10)
  expect_true(!is.null(res$weights_summary$ess))
  expect_true(res$weights_summary$ess <= res$n + 1e-8)
})
