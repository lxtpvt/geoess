ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) "")
if (!is.character(ofile) || length(ofile) != 1L || is.na(ofile)) ofile <- ""
base_dir <- if (nzchar(ofile)) dirname(ofile) else getwd()
pkg_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
force_source <- interactive() && nzchar(ofile)

if (!("geoess" %in% loadedNamespaces()) || force_source) {
  source(file.path(pkg_root, "R", "utils_checks.R"))
  source(file.path(pkg_root, "R", "utils_stability.R"))
  source(file.path(pkg_root, "R", "ess_core.R"))
  source(file.path(pkg_root, "R", "ess_misspec.R"))
  source(file.path(pkg_root, "R", "mixfit_mvnorm_em.R"))
  source(file.path(pkg_root, "R", "mix_mvnorm_calculus.R"))
  source(file.path(pkg_root, "R", "fisher_info.R"))
  source(file.path(pkg_root, "R", "stan_interface.R"))
  source(file.path(pkg_root, "R", "stan_geoess.R"))
  source(file.path(pkg_root, "R", "theta_ref.R"))
  source(file.path(pkg_root, "R", "hier.R"))
}

test_that("geoess_draws and theta_ref work with rstan", {
  testthat::skip_if_not_installed("rstan")

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1)

  stan_code <- "
  data {
    int<lower=0> N;
    vector[N] y;
    real<lower=0> sigma;
    real<lower=0> s0;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, s0);
    y ~ normal(mu, sigma);
  }
  generated quantities {
    vector[N] log_lik;
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu, sigma);
  }
  "

  y <- rnorm(10, mean = 0, sd = 2)
  data <- list(N = length(y), y = y, sigma = 2, s0 = 2)

  sm <- rstan::stan_model(model_code = stan_code)
  fit <- rstan::sampling(sm, data = data, iter = 200, chains = 1, refresh = 0, seed = 123)

  draws <- geoess_draws(fit, pars = "mu", as = "matrix")
  expect_true(is.matrix(draws))
  expect_true(nrow(draws) > 0)

  loglik <- geoess_loglik_draws(fit, var = "log_lik", as = "array")
  expect_true(is.array(loglik))

  ref <- geoess_theta_ref(fit, pars = "mu", ref = "posterior_mean")
  expect_true(is.numeric(ref))
  expect_equal(length(ref), 1)
})

test_that("geoess_draws and theta_ref work with cmdstanr", {
  testthat::skip_if_not_installed("cmdstanr")
  ver <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  if (is.null(ver)) testthat::skip("CmdStan not installed")

  stan_code <- "
  data {
    int<lower=0> N;
    vector[N] y;
    real<lower=0> sigma;
    real<lower=0> s0;
  }
  parameters {
    real mu;
  }
  model {
    mu ~ normal(0, s0);
    y ~ normal(mu, sigma);
  }
  generated quantities {
    vector[N] log_lik;
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu, sigma);
  }
  "

  y <- rnorm(10, mean = 0, sd = 2)
  data <- list(N = length(y), y = y, sigma = 2, s0 = 2)

  file <- cmdstanr::write_stan_file(stan_code)
  model <- cmdstanr::cmdstan_model(file)
  fit <- model$sample(data = data, iter_sampling = 200, iter_warmup = 100,
                      chains = 1, refresh = 0, seed = 123)

  draws <- geoess_draws(fit, pars = "mu", as = "matrix")
  expect_true(is.matrix(draws))
  expect_true(nrow(draws) > 0)

  loglik <- geoess_loglik_draws(fit, var = "log_lik", as = "array")
  expect_true(is.array(loglik))

  ref <- geoess_theta_ref(fit, pars = "mu", ref = "posterior_mean")
  expect_true(is.numeric(ref))
  expect_equal(length(ref), 1)
})
