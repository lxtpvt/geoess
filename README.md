# geoess

`geoess` is an R package for geometric effective sample size (ESS) in Bayesian
models. It treats ESS as a relative information quantity by comparing prior
curvature to the Fisher--Rao information contributed by one sampling unit.

The package supports:

- local ESS at user-chosen reference points,
- global ESS under prior, posterior, predictive, or design averaging,
- directional and spectral ESS for multivariate priors,
- hierarchical ESS for block and conditional borrowing, and
- sample-based workflows built from fitted mixtures or Stan draws.

The goal is to make prior informativeness and information borrowing measurable
in a way that stays coherent across conjugate, nonconjugate, multivariate, and
hierarchical settings.

## Installation

Install the development version from GitHub with:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("lxtpvt/geoess")
```

If you are working from a local checkout, you can also install with:

```r
install.packages(".", repos = NULL, type = "source")
```

## Key Functions

- `geoess_curvature()` for scalar/directional/eigen ESS
- `geoess_directional()` and `geoess_eigen()` for focused summaries
- `geoess_godambe()` for misspecified models (Godambe metric)
- `schur_theta_given_phi()` for hierarchical curvature (Schur complement)
- Normal-mean analytic module: `I1_normal_mu()`, `Ipi_normal_mu()`, `geoess_normal_mu()`
- Predictive ESS: `predictive_ess_normal_mu()` and `predictive_ess_discrete()`
- Stan integration: `geoess_draws()`, `geoess_theta_ref()`, `geoess_from_stan()`
- Mixture workflow: `geoess_mixfit()`, `geoess_automixfit()`, `geoess_Ipi_from_mix()`
- Fisher information: `geoess_I1_hessian()`, `geoess_I1_score()`
- Beta mixture tools: `mixbeta()`, `dmix()`, `pmix()`, `qmix()`, `rmix()`,
  `geoess_mixfit_beta()`, `geoess_automixfit_beta()`

## Examples

### Normal Mean with Known Variance

```r
library(geoess)

geoess_normal_mu(sigma = 2, s0 = 4)

demo_normal_mu(n = 50, mu_true = 1, sigma = 2, m0 = 0, s0 = 4, seed = 123)
```

### Directional and Eigen ESS

```r
I1 <- diag(c(2, 3))
Ipi <- diag(c(1, 6))

geoess_curvature(I1, Ipi)$ess
directional <- geoess_directional(I1, Ipi, v = c(2, -1))
values <- geoess_eigen(I1, Ipi)
list(directional = directional, eigen = values)
```

### Misspecification and Hierarchical Curvature

```r
I <- diag(c(2, 3))
J <- I
Ipi <- diag(c(1, 1))
geoess_godambe(I, J, Ipi)$ess

Htt <- matrix(2, 1, 1)
Htp <- matrix(1, 1, 1)
Hpp <- matrix(3, 1, 1)
schur_theta_given_phi(Htt, Htp, Hpp)
```

### Predictive ESS (Normal)

```r
predictive_ess_normal_mu(sigma = 2, s0 = 2, n0 = 1)
```

### Predictive ESS (Discrete Template)

```r
baseline <- c(0.5, 0.5)
target <- c(0.7, 0.3)
candidates <- rbind(c(0.55, 0.45), c(0.6, 0.4), c(0.7, 0.3))

predictive_ess_discrete(target, candidates, baseline, n_grid = 1:3)
```

### Visualization

```r
# 1D draws (hist + density)
draws <- cbind(rnorm(500))
geoess_plot_draws_1d(draws, method = "both")

# 2D draws (scatter or contour)
draws2 <- cbind(mu = rnorm(500), tau = rnorm(500, 1, 0.2))
geoess_plot_draws_2d(draws2, pars = c("mu", "tau"), method = "scatter")
geoess_plot_draws_2d(draws2, pars = c("mu", "tau"), method = "contour")

# Mixture fit visualization
mix <- geoess_automixfit(draws, family = "norm", K_grid = 1:3, criterion = "BIC")
geoess_plot_mix_1d(mix, dim_index = 1)
```

### Stan Integration (Sketch)

```r
# Requires rstan (imported) or cmdstanr (suggested) with CmdStan installed.
# The Hessian route needs loglik_fn; cmdstanr requires it explicitly.
# fit <- rstan::sampling(...) or cmdstanr::model$sample(...)
# prior_fit <- fit with N = 0 for prior draws, or use prior_draws matrix

loglik_fn <- function(theta, y, sigma) {
  sum(dnorm(y, mean = theta, sd = sigma, log = TRUE))
}

geoess_from_stan(
  fit,
  pars = "mu",
  n_obs = length(y),
  prior_draws_fit = prior_fit,
  theta_ref = "posterior_mean",
  I1_method = "hessian",
  loglik_fn = loglik_fn,
  y = y,
  sigma = sigma
)
```

### Global ESS (Expectation of Local ESS)

```r
sigma <- 2
s0 <- 3
I1 <- I1_normal_mu(sigma)
Ipi <- Ipi_normal_mu(s0)

set.seed(1)
draws <- matrix(rnorm(1000, mean = 0, sd = s0), ncol = 1)

# Prior-averaged ESS under Q = prior
geoess_global(draws, I1 = I1, Ipi = Ipi, method = "trace")
```

### Beta Mixture Example

```r
mix <- mixbeta(c(0.6, 2, 5), c(0.4, 5, 2), param = "ab")
dmix(mix, c(0.2, 0.5))
pmix(mix, 0.5)
qmix(mix, 0.5)
rmix(mix, 5)
```

### Beta Mixture Fit (EM)

```r
set.seed(123)
draws <- rbeta(300, 2, 6)
fit <- geoess_mixfit_beta(draws, K = 1)
attr(fit, "logLik")
```
