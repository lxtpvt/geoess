#' Reference Point Selection for Local ESS
#'
#' ESS is local and depends on the evaluation point \eqn{\theta^*}. This function
#' selects \eqn{\theta^*} (or \eqn{p^*} for probability parameters) using
#' prior/posterior mean, median, MAP, pseudo-true (KL projection), or
#' target-matching rules.
#'
#' @param source A Stan fit, draws matrix, \code{geoess_mix} object, or
#'   \code{geoess_spec} list containing \code{prior_source} and
#'   \code{posterior_source} entries.
#' @param pars Optional character vector of parameter names to extract from draws.
#' @param ref Character reference choice.
#' @param domain Character. "real" or "prob" (probabilities in (0,1)).
#' @param target Numeric scalar target ESS (for \code{ref = "target_match"}).
#' @param target_method Character ESS scalarization method.
#' @param v Optional direction vector for \code{target_method = "dir"}.
#' @param lower Optional lower bounds (or lower probability).
#' @param upper Optional upper bounds (or upper probability).
#' @param init Optional initial guess.
#' @param data Optional data for log-likelihood or pseudo-true objectives.
#' @param model Optional list with model-specific functions/metadata.
#' @param loglik_fn Optional log-likelihood function.
#' @param logprior_fn Optional log-prior function.
#' @param mix Optional mixture object for prior/posterior density.
#' @param optimizer Character. "nlminb" or "optim".
#' @param control Optional optimizer control list.
#' @param ... Passed to \code{loglik_fn} / \code{logprior_fn}.
#'
#' @return Named numeric reference vector with attributes:
#' \code{ref}, \code{objective}, \code{converged}, and \code{details}.
#' @export
geoess_theta_ref <- function(source,
                             pars = NULL,
                             ref = c("prior_mode", "posterior_mode",
                                     "prior_mean", "posterior_mean",
                                     "prior_median", "posterior_median",
                                     "pseudo_true", "target_match"),
                             domain = c("real", "prob"),
                             target = NULL,
                             target_method = c("trace", "dir", "eig_min", "eig_mean"),
                             v = NULL,
                             lower = NULL, upper = NULL,
                             init = NULL,
                             data = NULL,
                             model = NULL,
                             loglik_fn = NULL,
                             logprior_fn = NULL,
                             mix = NULL,
                             optimizer = c("nlminb", "optim"),
                             control = list(),
                             ...) {
  ref <- match.arg(ref)
  domain <- match.arg(domain)
  target_method <- match.arg(target_method)
  optimizer <- match.arg(optimizer)

  ref_type <- if (grepl("^prior_", ref)) "prior" else if (grepl("^posterior_", ref)) "posterior" else NULL
  details <- list()

  if (inherits(source, "geoess_spec")) {
    spec <- source
    if (!is.list(spec)) stop("geoess_spec must be a list")
    if (is.null(data)) data <- spec$data
    if (is.null(model)) model <- spec$model
    if (is.null(loglik_fn)) loglik_fn <- spec$loglik_fn
    if (is.null(logprior_fn)) logprior_fn <- spec$logprior_fn
    if (is.null(mix)) mix <- spec$mix

    if (ref_type == "prior") {
      source <- spec$prior_source
    } else if (ref_type == "posterior") {
      source <- spec$posterior_source
    } else {
      source <- spec$posterior_source
    }

    if (is.null(source)) source <- spec$prior_source
    if (is.null(source)) source <- spec$source
    if (is.null(source)) {
      stop("geoess_spec does not contain a source for ref = ", ref)
    }
  }

  if (ref %in% c("prior_mean", "prior_median", "posterior_mean", "posterior_median")) {
    mix_use <- mix
    if (ref_type == "posterior" &&
        !inherits(source, "geoess_mix") &&
        !inherits(source, "betaMix")) {
      mix_use <- NULL
    }
    res <- .geoess_ref_mean_median(
      source = source,
      mix = mix_use,
      pars = pars,
      type = ref_type,
      median = grepl("median$", ref),
      domain = domain,
      lower = lower,
      upper = upper
    )
    out <- res$theta
    details <- c(details, res$details)
    attr(out, "ref") <- ref
    attr(out, "objective") <- NA_real_
    attr(out, "converged") <- TRUE
    attr(out, "details") <- details
    return(out)
  }

  if (ref == "prior_mode") {
    res <- .geoess_ref_mode(
      source = source,
      mix = mix,
      pars = pars,
      domain = domain,
      lower = lower,
      upper = upper,
      init = init,
      optimizer = optimizer,
      control = control,
      logprior_fn = logprior_fn,
      ...
    )
    out <- res$theta
    details <- c(details, res$details)
    attr(out, "ref") <- ref
    attr(out, "objective") <- res$objective
    attr(out, "converged") <- res$converged
    attr(out, "details") <- details
    return(out)
  }

  if (ref == "posterior_mode") {
    res <- .geoess_ref_posterior_mode(
      source = source,
      mix = mix,
      pars = pars,
      domain = domain,
      lower = lower,
      upper = upper,
      init = init,
      optimizer = optimizer,
      control = control,
      loglik_fn = loglik_fn,
      logprior_fn = logprior_fn,
      ...
    )
    out <- res$theta
    details <- c(details, res$details)
    attr(out, "ref") <- ref
    attr(out, "objective") <- res$objective
    attr(out, "converged") <- res$converged
    attr(out, "details") <- details
    return(out)
  }

  if (ref == "pseudo_true") {
    res <- .geoess_ref_pseudo_true(
      source = source,
      pars = pars,
      domain = domain,
      lower = lower,
      upper = upper,
      init = init,
      optimizer = optimizer,
      control = control,
      loglik_fn = loglik_fn,
      data = data,
      model = model,
      ...
    )
    out <- res$theta
    details <- c(details, res$details)
    attr(out, "ref") <- ref
    attr(out, "objective") <- res$objective
    attr(out, "converged") <- res$converged
    attr(out, "details") <- details
    return(out)
  }

  res <- .geoess_ref_target_match(
    source = source,
    mix = mix,
    pars = pars,
    domain = domain,
    target = target,
    target_method = target_method,
    v = v,
    lower = lower,
    upper = upper,
    init = init,
    optimizer = optimizer,
    control = control,
    loglik_fn = loglik_fn,
    logprior_fn = logprior_fn,
    data = data,
    model = model,
    ...
  )
  out <- res$theta
  details <- c(details, res$details)
  attr(out, "ref") <- ref
  attr(out, "objective") <- res$objective
  attr(out, "converged") <- res$converged
  attr(out, "details") <- details
  out
}

.geoess_ref_mean_median <- function(source, mix, pars, type, median,
                                    domain, lower, upper) {
  details <- list(type = type)

  if (inherits(source, "geoess_mix") || inherits(source, "betaMix")) {
    mix <- source
  }

  if (!is.null(mix)) {
    if (inherits(mix, "betaMix")) {
      if (median) {
        theta <- qmix(mix, 0.5)
      } else {
        w <- mix["w", ]
        a <- mix["a", ]
        b <- mix["b", ]
        theta <- sum(w * (a / (a + b)))
      }
      return(list(theta = as.numeric(theta), details = details))
    }

    if (inherits(mix, "geoess_mix")) {
      w <- as.numeric(mix$weights)
      w <- w / sum(w)
      mu <- do.call(rbind, mix$mu)
      if (!median) {
        theta <- as.numeric(crossprod(w, mu))
        return(list(theta = theta, details = details))
      }

      if (ncol(mu) != 1L) {
        stop("Median for multivariate mixtures requires draws")
      }

      sd <- vapply(mix$Sigma, function(S) sqrt(S[1, 1]), numeric(1))
      cdf_mix <- function(x) sum(w * stats::pnorm(x, mean = mu[, 1], sd = sd))
      rng <- .geoess_bracket_mix(mu[, 1], sd, lower = lower, upper = upper)
      theta <- stats::uniroot(function(x) cdf_mix(x) - 0.5, rng)$root
      return(list(theta = as.numeric(theta), details = details))
    }
  }

  if (is.null(type)) stop("prior/posterior type is required for draws")
  draws <- .geoess_as_draws_matrix(source, pars = pars, type = type)
  if (median) {
    theta <- apply(draws, 2, stats::median)
  } else {
    theta <- colMeans(draws)
  }
  list(theta = as.numeric(theta), details = details)
}

.geoess_ref_mode <- function(source, mix, pars, domain, lower, upper, init,
                             optimizer, control, logprior_fn, ...) {
  details <- list()
  if (inherits(source, "geoess_mix") || inherits(source, "betaMix")) {
    mix <- source
  }

  if (is.null(logprior_fn)) {
    if (inherits(mix, "geoess_mix")) {
      logprior_fn <- function(theta) geoess_mix_logdens(mix, theta)
    } else if (inherits(mix, "betaMix")) {
      logprior_fn <- function(theta) dmix(mix, theta, log = TRUE)
    } else {
      stop("logprior_fn or mix is required for prior_mode")
    }
  }

  if (is.null(init)) {
    init <- .geoess_default_init(source = source, mix = mix, pars = pars, domain = domain)
  }

  obj <- function(theta) -logprior_fn(theta, ...)
  res <- .geoess_optimize(obj, init, domain, lower, upper, optimizer, control)
  details$optimizer <- res$optimizer
  list(theta = res$par, objective = res$objective, converged = res$converged, details = details)
}

.geoess_ref_posterior_mode <- function(source, mix, pars, domain, lower, upper,
                                       init, optimizer, control,
                                       loglik_fn, logprior_fn, ...) {
  details <- list()
  if (is.null(loglik_fn)) {
    map_val <- .geoess_theta_ref_map(source, pars, ...)
    if (!is.null(map_val)) {
      details$optimizer <- "stan"
      return(list(theta = map_val, objective = NA_real_, converged = TRUE, details = details))
    }
    stop("loglik_fn is required for posterior_mode")
  }

  if (inherits(source, "geoess_mix") || inherits(source, "betaMix")) {
    mix <- source
  }

  if (is.null(logprior_fn)) {
    if (inherits(mix, "geoess_mix")) {
      logprior_fn <- function(theta) geoess_mix_logdens(mix, theta)
    } else if (inherits(mix, "betaMix")) {
      logprior_fn <- function(theta) dmix(mix, theta, log = TRUE)
    } else {
      logprior_fn <- function(theta) 0
    }
  }

  if (is.null(init)) {
    init <- .geoess_default_init(source = source, mix = mix, pars = pars, domain = domain)
  }

  obj <- function(theta) -(loglik_fn(theta, ...) + logprior_fn(theta, ...))
  res <- .geoess_optimize(obj, init, domain, lower, upper, optimizer, control)
  details$optimizer <- res$optimizer
  list(theta = res$par, objective = res$objective, converged = res$converged, details = details)
}

.geoess_ref_pseudo_true <- function(source, pars, domain, lower, upper, init,
                                    optimizer, control, loglik_fn, data, model, ...) {
  if (is.null(loglik_fn) && (is.null(model) || is.null(model$loglik_point))) {
    stop("loglik_fn or model$loglik_point is required for pseudo_true")
  }

  mode <- if (!is.null(model) && !is.null(model$pseudo_true_mode)) {
    model$pseudo_true_mode
  } else if (!is.null(model) && !is.null(model$mode)) {
    model$mode
  } else {
    "empirical"
  }

  details <- list(pseudo_true_mode = mode)

  if (is.null(init)) {
    init <- .geoess_default_init(source = source, pars = pars, domain = domain)
  }

  if (mode == "mc") {
    if (is.null(model$r_true) || is.null(model$loglik_point)) {
      stop("model$r_true and model$loglik_point are required for mc pseudo_true")
    }
    n_mc <- if (!is.null(model$n_mc)) model$n_mc else if (!is.null(data) && !is.null(data$n_mc)) data$n_mc else 1000
    x <- model$r_true(n_mc)

    obj <- function(theta) {
      ll <- model$loglik_point(theta, x)
      -mean(ll)
    }
  } else {
    obj <- function(theta) -loglik_fn(theta, ...)
  }

  res <- .geoess_optimize(obj, init, domain, lower, upper, optimizer, control)
  list(theta = res$par, objective = res$objective, converged = res$converged, details = details)
}

.geoess_ref_target_match <- function(source, mix, pars, domain, target, target_method,
                                     v, lower, upper, init, optimizer, control,
                                     loglik_fn, logprior_fn, data, model, ...) {
  if (is.null(target) || length(target) != 1L || !is.finite(target)) {
    stop("target must be a finite numeric scalar")
  }

  if (is.null(model)) model <- list()
  psd <- if (!is.null(model$psd)) model$psd else "none"

  if (inherits(source, "geoess_mix") || inherits(source, "betaMix")) {
    mix <- source
  }

  I1_fn <- if (!is.null(model$I1_fn)) model$I1_fn else NULL

  prob_bounds <- if (domain == "prob") .geoess_bracket_prob(lower, upper) else NULL

  ess_scalar <- function(theta) {
    I1 <- if (!is.null(I1_fn)) {
      geoess_I1_from_fn(theta, I1_fn)
    } else if (!is.null(loglik_fn)) {
      n_obs <- .geoess_n_obs(data, model)
      geoess_I1_from_hessian(theta, loglik_fn, n_obs, ...)
    } else {
      stop("I1_fn or loglik_fn is required for target_match")
    }

    Ipi <- NULL
    if (inherits(mix, "geoess_mix")) {
      Ipi <- geoess_Ipi_from_mix(mix, theta, psd = psd)
    } else {
      if (is.null(logprior_fn)) {
        if (inherits(mix, "betaMix")) {
          logprior_fn <- function(x) dmix(mix, x, log = TRUE)
        } else {
          stop("logprior_fn or mix is required for target_match")
        }
      }
      if (domain == "prob" && length(theta) == 1L) {
        Ipi <- .geoess_hess_prob_scalar(logprior_fn, theta, prob_bounds[1], prob_bounds[2], ...)
      } else {
        Ipi <- -numDeriv::hessian(logprior_fn, theta, ...)
        Ipi <- .geoess_symmetrize(Ipi)
      }
    }

    if (any(!is.finite(as.numeric(I1))) || any(!is.finite(as.numeric(Ipi)))) return(NA_real_)

    if (target_method == "trace") {
      geoess_ess_trace(I1, Ipi)
    } else if (target_method == "dir") {
      if (is.null(v)) stop("v is required for target_method = 'dir'")
      geoess_ess_dir(I1, Ipi, v)
    } else {
      eigs <- geoess_ess_eigs(I1, Ipi)
      if (target_method == "eig_min") min(eigs) else mean(eigs)
    }
  }

  details <- list(target = target, target_method = target_method, psd = psd)

  if (domain == "prob") {
    rng <- prob_bounds
    eps <- min(1e-4, 0.1 * (rng[2] - rng[1]))
    grid <- seq(rng[1] + eps, rng[2] - eps, length.out = 200)
    vals <- vapply(grid, function(p) ess_scalar(p) - target, numeric(1))
    idx <- integer(0)
    for (i in seq_len(length(vals) - 1L)) {
      if (is.finite(vals[i]) && is.finite(vals[i + 1L]) && vals[i] * vals[i + 1L] <= 0) {
        idx <- i
        break
      }
    }
    if (length(idx) == 0) {
      warning("Unable to bracket target; returning NA")
      return(list(theta = NA_real_, objective = NA_real_, converged = FALSE, details = details))
    }
    lower_b <- grid[idx[1]]
    upper_b <- grid[idx[1] + 1]
    sol <- stats::uniroot(function(p) ess_scalar(p) - target,
                          lower = lower_b, upper = upper_b, tol = 1e-10)
    return(list(theta = sol$root, objective = sol$f.root, converged = TRUE, details = details))
  }

  path_fn <- if (!is.null(model$path_fn)) model$path_fn else NULL
  if (!is.null(path_fn)) {
    t_bounds <- if (!is.null(model$t_bounds)) model$t_bounds else c(0, 1)
    sol <- stats::uniroot(function(t) ess_scalar(path_fn(t)) - target,
                          lower = t_bounds[1], upper = t_bounds[2])
    return(list(theta = path_fn(sol$root), objective = sol$f.root, converged = TRUE, details = details))
  }

  if (is.null(init)) stop("init is required for target_match in multi-dim")
  obj <- function(theta) {
    val <- ess_scalar(theta) - target
    val^2
  }
  res <- .geoess_optimize(obj, init, domain, lower, upper, optimizer, control)
  list(theta = res$par, objective = res$objective, converged = res$converged, details = details)
}

.geoess_optimize <- function(fn, init, domain, lower, upper, optimizer, control) {
  if (domain == "prob") {
    low <- if (is.null(lower)) 0 else lower
    high <- if (is.null(upper)) 1 else upper
    if (low <= 0 || high >= 1 || low >= high) stop("Invalid bounds for prob domain")
    trans <- function(u) low + (high - low) * .geoess_invlogit(u)
    invtrans <- function(p) qlogis((p - low) / (high - low))
    u0 <- if (is.null(init)) 0 else invtrans(init)
    obj_u <- function(u) fn(trans(u))
    if (optimizer == "nlminb") {
      fit <- do.call(stats::nlminb, c(list(start = u0, objective = obj_u), control))
      return(list(par = trans(fit$par), objective = fit$objective,
                  converged = fit$convergence == 0, optimizer = "nlminb"))
    }
    fit <- stats::optim(u0, obj_u, method = "BFGS", control = control)
    return(list(par = trans(fit$par), objective = fit$value,
                converged = fit$convergence == 0, optimizer = "optim"))
  }

  init <- as.numeric(init)
  if (optimizer == "nlminb") {
    args <- list(start = init, objective = fn, control = control)
    if (!is.null(lower)) args$lower <- lower
    if (!is.null(upper)) args$upper <- upper
    fit <- do.call(stats::nlminb, args)
    return(list(par = fit$par, objective = fit$objective,
                converged = fit$convergence == 0, optimizer = "nlminb"))
  }

  method <- if (!is.null(lower) || !is.null(upper)) "L-BFGS-B" else "BFGS"
  fit <- stats::optim(init, fn, method = method, lower = lower, upper = upper, control = control)
  list(par = fit$par, objective = fit$value,
       converged = fit$convergence == 0, optimizer = "optim")
}

.geoess_default_init <- function(source, mix = NULL, pars = NULL, domain = "real") {
  if (domain == "prob") return(0.5)
  if (!is.null(mix)) {
    if (inherits(mix, "geoess_mix")) {
      w <- as.numeric(mix$weights)
      w <- w / sum(w)
      mu <- do.call(rbind, mix$mu)
      return(as.numeric(crossprod(w, mu)))
    }
    if (inherits(mix, "betaMix")) {
      w <- mix["w", ]
      a <- mix["a", ]
      b <- mix["b", ]
      return(sum(w * (a / (a + b))))
    }
  }
  if (!is.null(source)) {
    draws <- tryCatch(.geoess_as_draws_matrix(source, pars = pars, type = "posterior"),
                      error = function(e) NULL)
    if (!is.null(draws)) return(colMeans(draws))
  }
  0
}

.geoess_as_draws_matrix <- function(source, pars, type = c("prior", "posterior")) {
  type <- match.arg(type)
  if (inherits(source, "geoess_spec")) {
    src <- if (type == "prior") source$prior_source else source$posterior_source
    if (is.null(src)) stop("Missing ", type, " source in geoess_spec")
    source <- src
  }
  if (inherits(source, "geoess_mix") || inherits(source, "betaMix")) {
    stop("Mixture object does not contain draws")
  }
  geoess_draws(source, pars = pars, type = type, as = "matrix")
}

.geoess_bracket_mix <- function(mu, sd, lower, upper) {
  if (!is.null(lower) && !is.null(upper)) return(c(lower, upper))
  lo <- min(mu - 6 * sd)
  hi <- max(mu + 6 * sd)
  c(lo, hi)
}

.geoess_bracket_prob <- function(lower, upper) {
  lo <- if (is.null(lower)) 1e-4 else lower
  hi <- if (is.null(upper)) 1 - 1e-4 else upper
  if (lo <= 0 || hi >= 1 || lo >= hi) stop("Invalid bounds for prob domain")
  c(lo, hi)
}

.geoess_hess_prob_scalar <- function(logprior_fn, p, lower, upper, ...) {
  p <- as.numeric(p)
  if (!is.finite(p) || p <= lower || p >= upper) return(matrix(NA_real_, nrow = 1L))

  margin <- min(p - lower, upper - p)
  if (!is.finite(margin) || margin <= 0) return(matrix(NA_real_, nrow = 1L))

  h <- min(1e-4, 0.49 * margin)
  if (!is.finite(h) || h <= 0) return(matrix(NA_real_, nrow = 1L))

  f0 <- logprior_fn(p, ...)
  f1 <- logprior_fn(p + h, ...)
  f2 <- logprior_fn(p - h, ...)
  if (!is.finite(f0) || !is.finite(f1) || !is.finite(f2)) {
    return(matrix(NA_real_, nrow = 1L))
  }

  hess <- (f1 - 2 * f0 + f2) / (h^2)
  matrix(-hess, nrow = 1L)
}

.geoess_n_obs <- function(data, model) {
  if (!is.null(model) && !is.null(model$n_obs)) return(model$n_obs)
  if (!is.null(data) && !is.null(data$n_obs)) return(data$n_obs)
  if (is.list(data) && !is.null(data$y)) return(length(data$y))
  if (is.numeric(data)) return(length(data))
  stop("n_obs could not be inferred; provide model$n_obs or data$n_obs")
}

.geoess_invlogit <- function(x) {
  1 / (1 + exp(-x))
}
