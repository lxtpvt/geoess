ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "mix_mvnorm_calculus.R"))
  source(file.path(pkg_root, "R", "fisher_info.R"))
  source(file.path(pkg_root, "R", "stan_interface.R"))
  source(file.path(pkg_root, "R", "theta_ref.R"))
}

test_that("geoess_theta_ref returns mean/median from draws", {
  set.seed(101)
  draws <- matrix(rnorm(1000, mean = 2, sd = 1), ncol = 1)
  colnames(draws) <- "theta"

  ref_mean <- geoess_theta_ref(draws, pars = "theta", ref = "prior_mean")
  ref_med <- geoess_theta_ref(draws, pars = "theta", ref = "prior_median")

  expect_equal(as.numeric(ref_mean), mean(draws), tolerance = 1e-12)
  expect_equal(as.numeric(ref_med), stats::median(draws), tolerance = 1e-12)
})

test_that("geoess_theta_ref prior_mode matches MVN component mean", {
  mix <- structure(
    list(
      family = "mvnorm",
      weights = 1,
      mu = list(c(1.5, -0.5)),
      Sigma = list(diag(c(2, 3)))
    ),
    class = "geoess_mix"
  )

  ref <- geoess_theta_ref(mix, ref = "prior_mode")
  expect_equal(as.numeric(ref), c(1.5, -0.5), tolerance = 1e-6)
})

test_that("geoess_theta_ref posterior_mode matches Normal MAP", {
  set.seed(123)
  y <- rnorm(20, mean = 1, sd = 2)
  sigma <- 2
  m0 <- -0.5
  s0 <- 1.5

  loglik_fn <- function(theta, y, sigma, ...) {
    sum(dnorm(y, mean = theta, sd = sigma, log = TRUE))
  }
  logprior_fn <- function(theta, m0, s0, ...) {
    dnorm(theta, mean = m0, sd = s0, log = TRUE)
  }

  map_true <- (length(y) * mean(y) / sigma^2 + m0 / s0^2) /
    (length(y) / sigma^2 + 1 / s0^2)

  ref <- geoess_theta_ref(
    source = matrix(y, ncol = 1),
    ref = "posterior_mode",
    init = 0,
    loglik_fn = loglik_fn,
    logprior_fn = logprior_fn,
    y = y,
    sigma = sigma,
    m0 = m0,
    s0 = s0
  )

  expect_equal(as.numeric(ref), map_true, tolerance = 1e-5)
})

test_that("geoess_theta_ref pseudo_true matches empirical MLE", {
  set.seed(321)
  y <- rnorm(30, mean = 0.3, sd = 1.5)
  sigma <- 1.5

  loglik_fn <- function(theta, y, sigma) {
    sum(dnorm(y, mean = theta, sd = sigma, log = TRUE))
  }

  ref <- geoess_theta_ref(
    source = matrix(y, ncol = 1),
    ref = "pseudo_true",
    init = 0,
    loglik_fn = loglik_fn,
    y = y,
    sigma = sigma
  )

  expect_equal(as.numeric(ref), mean(y), tolerance = 1e-5)
})

test_that("geoess_theta_ref target_match solves Beta ESS", {
  a <- 5
  b <- 4
  target <- 8

  logprior_fn <- function(p, a, b) dbeta(p, a, b, log = TRUE)
  I1_fn <- function(p) matrix(1 / (p * (1 - p)), nrow = 1)

  ref <- geoess_theta_ref(
    source = matrix(0, ncol = 1),
    ref = "target_match",
    domain = "prob",
    target = target,
    lower = 0.01,
    upper = 0.99,
    model = list(I1_fn = I1_fn),
    logprior_fn = logprior_fn,
    a = a,
    b = b
  )

  ref_val <- as.numeric(ref)
  I1 <- I1_fn(ref_val)
  Ipi <- .geoess_hess_prob_scalar(function(p) logprior_fn(p, a, b),
                                  ref_val, lower = 0.01, upper = 0.99)
  ess <- geoess_trace(I1, Ipi)

  expect_true(is.finite(ref_val))
  expect_equal(as.numeric(ess), target, tolerance = 1e-5)
})
