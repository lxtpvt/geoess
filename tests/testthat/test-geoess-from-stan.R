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

test_that("geoess_from_stan works with rstan", {
  testthat::skip_if_not_installed("rstan")

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 1)
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

  set.seed(1)
  y <- rnorm(20, mean = 0, sd = 2)
  data <- list(N = length(y), y = y, sigma = 2, s0 = 4)
  prior_data <- list(N = 0, y = numeric(0), sigma = 2, s0 = 4)

  sm <- rstan::stan_model(model_code = stan_code)
  fit <- rstan::sampling(sm, data = data, iter = 200, chains = 1, refresh = 0, seed = 123)
  prior_fit <- rstan::sampling(sm, data = prior_data, iter = 200, chains = 1, refresh = 0, seed = 124)

  loglik_fn <- function(theta, y, sigma) {
    sum(dnorm(y, mean = theta, sd = sigma, log = TRUE))
  }

  set.seed(123)
  res <- geoess_from_stan(
    fit,
    pars = "mu",
    n_obs = length(y),
    prior_draws_fit = prior_fit,
    mixture_family = "norm",
    mixture_K_grid = 1:2,
    mixture_criterion = "BIC",
    theta_ref = "posterior_mean",
    I1_method = "hessian",
    loglik_fn = loglik_fn,
    method = "trace",
    y = y,
    sigma = 2
  )

  expect_true(is.numeric(res$ess))
  expect_true(abs(res$ess - (data$sigma^2 / data$s0^2)) < 0.2)
})

test_that("geoess_from_stan works with cmdstanr", {
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

  set.seed(1)
  y <- rnorm(20, mean = 0, sd = 2)
  data <- list(N = length(y), y = y, sigma = 2, s0 = 4)
  prior_data <- list(N = 0, y = numeric(0), sigma = 2, s0 = 4)

  file <- cmdstanr::write_stan_file(stan_code)
  model <- cmdstanr::cmdstan_model(file)
  fit <- model$sample(data = data, iter_sampling = 200, iter_warmup = 100,
                      chains = 1, refresh = 0, seed = 123)
  prior_fit <- model$sample(data = prior_data, iter_sampling = 200, iter_warmup = 100,
                            chains = 1, refresh = 0, seed = 124)

  loglik_fn <- function(theta, y, sigma) {
    sum(dnorm(y, mean = theta, sd = sigma, log = TRUE))
  }

  set.seed(123)
  res <- geoess_from_stan(
    fit,
    pars = "mu",
    n_obs = length(y),
    prior_draws_fit = prior_fit,
    mixture_family = "norm",
    mixture_K_grid = 1:2,
    mixture_criterion = "BIC",
    theta_ref = "posterior_mean",
    I1_method = "hessian",
    loglik_fn = loglik_fn,
    method = "trace",
    y = y,
    sigma = 2
  )

  expect_true(is.numeric(res$ess))
  expect_true(abs(res$ess - (data$sigma^2 / data$s0^2)) < 0.2)
})
