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
  source(file.path(pkg_root, "R", "mix_mvnorm_calculus.R"))
  source(file.path(pkg_root, "R", "local_mixture_ess.R"))
  source(file.path(pkg_root, "R", "global_ess.R"))
  source(file.path(pkg_root, "R", "global_ess_result.R"))
}

test_that("geoess_global prior-averaged ESS matches constant", {
  sigma <- 2
  s0 <- 3
  I1 <- I1_normal_mu(sigma)
  Ipi <- Ipi_normal_mu(s0)

  set.seed(1)
  draws <- matrix(rnorm(1000, mean = 0, sd = s0), ncol = 1)

  res <- geoess_global(draws, I1 = I1, Ipi = Ipi, method = "trace")
  expect_equal(res$estimate, sigma^2 / s0^2, tolerance = 1e-10)
})

test_that("geoess_global posterior-averaged ESS matches constant", {
  sigma <- 2
  s0 <- 3
  I1 <- I1_normal_mu(sigma)
  Ipi <- Ipi_normal_mu(s0)

  n <- 20
  y <- rnorm(n, mean = 1, sd = sigma)
  v_post <- 1 / (n / sigma^2 + 1 / s0^2)
  m_post <- v_post * (n * mean(y) / sigma^2)

  set.seed(2)
  draws <- matrix(rnorm(800, mean = m_post, sd = sqrt(v_post)), ncol = 1)

  res <- geoess_global(draws, I1 = I1, Ipi = Ipi, method = "trace")
  expect_equal(res$estimate, sigma^2 / s0^2, tolerance = 1e-10)
})

test_that("geoess_global design-based ESS works", {
  sigma <- 2
  s0 <- 3
  I1 <- I1_normal_mu(sigma)
  Ipi <- Ipi_normal_mu(s0)

  rQ <- function(n) matrix(rnorm(n, mean = 0, sd = s0), ncol = 1)
  res <- geoess_global_design(rQ, I1 = I1, Ipi = Ipi, nsim = 500)

  expect_equal(res$estimate, sigma^2 / s0^2, tolerance = 1e-10)
})

test_that("geoess_global uses weighted hard-assignment ESS for mixtures by default", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.7, 0.3),
      mu = list(c(-4), c(4)),
      Sigma = list(matrix(1, 1, 1), matrix(0.25, 1, 1)),
      prec = list(matrix(1, 1, 1), matrix(4, 1, 1))
    ),
    class = "geoess_mix"
  )
  draws <- matrix(c(-4, 4), ncol = 1)
  I1 <- matrix(1, 1, 1)

  res <- geoess_global(draws, I1 = I1, Ipi = mix, method = "trace",
                       return_samples = TRUE)

  expect_equal(res$ess_samples, c(0.7, 1.2), tolerance = 1e-10)
  expect_equal(res$estimate, mean(c(0.7, 1.2)), tolerance = 1e-10)
  expect_equal(res$meta$mixture_attribution, "hard_resp")
})

test_that("geoess_global PSD handling for full-mixture curvature remains available", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = c(0.5, 0.5),
      mu = list(c(-3), c(3)),
      Sigma = list(matrix(0.2^2, 1, 1), matrix(0.2^2, 1, 1))
    ),
    class = "geoess_mix"
  )
  draws <- matrix(seq(-1, 1, length.out = 50), ncol = 1)
  I1 <- matrix(1, 1, 1)

  expect_warning(
    res_none <- geoess_global(draws, I1 = I1, Ipi = mix, method = "trace",
                              mixture_attribution = "mix",
                              psd = "none", return_samples = TRUE)
  )
  expect_true(any(res_none$ess_samples < 0))

  res_clip <- geoess_global(draws, I1 = I1, Ipi = mix, method = "trace",
                            mixture_attribution = "mix",
                            psd = "clip", return_samples = TRUE)
  expect_true(all(res_clip$ess_samples >= -1e-10))
})
