#' Geometric ESS from Stan Fits
#'
#' Fits a mixture to prior draws, computes curvature and Fisher information,
#' and returns geometric ESS.
#'
#' @param fit Stan fit object with posterior draws.
#' @param pars Character vector of parameters.
#' @param n_obs Integer number of observations.
#' @param prior_draws_fit Optional Stan fit providing prior draws.
#' @param prior_draws Optional matrix/array of prior draws.
#' @param mixture_family Character. "mvnorm" or "norm".
#' @param mixture_K_grid Integer vector of candidate K values.
#' @param mixture_criterion Character. "BIC" or "AIC".
#' @param theta_ref Character or numeric. Reference point choice (e.g.,
#'   "posterior_mean", "posterior_mode", "prior_mean", "prior_mode",
#'   "posterior_median", "prior_median", "pseudo_true", "target_match") or
#'   a numeric vector to use directly.
#' @param theta_ref_args Optional list of arguments passed to
#'   \code{geoess_theta_ref} (e.g., \code{target}, \code{domain}, \code{init}).
#' @param I1_method Character. "hessian" or "score".
#' @param loglik_fn Optional log-likelihood function for Hessian route.
#' @param score_fn Optional score function for score route.
#' @param method Character. ESS method "trace", "eig", or "dir".
#' @param v Optional direction for directional ESS.
#' @param ... Passed to \code{loglik_fn} or \code{score_fn}.
#'
#' @details
#' The Hessian route requires \code{loglik_fn}. For rstan, a fallback to
#' \code{rstan::log_prob} is attempted only when the full parameter vector
#' is supplied; cmdstanr requires \code{loglik_fn}. For Stan fits, \code{pars}
#' should be parameter names. Reference points are selected with
#' \code{geoess_theta_ref}. For mixture priors, component-wise ESS values and a
#' weighted ESS are returned in \code{ess_components} and \code{ess_weighted}
#' when \code{method = "trace"}.
#'
#' @return A list with ESS outputs, curvature matrices, and metadata.
#' @export
geoess_from_stan <- function(fit,
                             pars,
                             n_obs,
                             prior_draws_fit = NULL,
                             prior_draws = NULL,
                             mixture_family = c("mvnorm", "norm"),
                             mixture_K_grid = 1:5,
                             mixture_criterion = c("BIC", "AIC"),
                             theta_ref = c("posterior_mean", "posterior_mode",
                                           "prior_mean", "prior_mode",
                                           "posterior_median", "prior_median",
                                           "pseudo_true", "target_match"),
                             theta_ref_args = list(),
                             I1_method = c("hessian", "score"),
                             loglik_fn = NULL,
                             score_fn = NULL,
                             method = c("trace", "eig", "dir"),
                             v = NULL,
                             ...) {
  if (missing(pars)) stop("pars must be provided")
  if (is.numeric(pars) &&
      (inherits(fit, "stanfit") || inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanVB"))) {
    stop("pars must be character names for Stan fits")
  }
  if (!is.numeric(pars)) pars <- as.character(pars)
  if (!is.numeric(n_obs) || length(n_obs) != 1L || n_obs <= 0) {
    stop("n_obs must be positive")
  }

  mixture_family <- match.arg(mixture_family)
  mixture_criterion <- match.arg(mixture_criterion)
  I1_method <- match.arg(I1_method)
  method <- match.arg(method)
  if (is.null(theta_ref_args)) theta_ref_args <- list()
  if (!is.list(theta_ref_args)) stop("theta_ref_args must be a list")

  if (is.null(prior_draws) && is.null(prior_draws_fit)) {
    stop("Provide prior_draws or prior_draws_fit")
  }

  if (!is.null(prior_draws)) {
    prior_mat <- geoess_draws(prior_draws, pars = pars, as = "matrix")
  } else {
    prior_mat <- geoess_draws(prior_draws_fit, pars = pars, type = "prior", as = "matrix")
  }

  mix <- geoess_automixfit(
    prior_mat,
    family = mixture_family,
    K_grid = mixture_K_grid,
    criterion = mixture_criterion
  )

  if (is.numeric(theta_ref)) {
    theta_hat <- as.numeric(theta_ref)
    if (length(theta_hat) != length(pars)) {
      stop("theta_ref length does not match pars")
    }
    if (is.numeric(pars) && !is.null(colnames(prior_mat))) {
      names(theta_hat) <- colnames(prior_mat)
    } else {
      names(theta_hat) <- pars
    }
  } else {
    if (length(theta_ref) == 1L) {
      if (theta_ref == "mean") theta_ref <- "posterior_mean"
      if (theta_ref == "map") theta_ref <- "posterior_mode"
    }
    theta_ref <- match.arg(theta_ref)
    if (grepl("^prior_", theta_ref)) {
      source_ref <- prior_mat
      pars_ref <- NULL
    } else {
      source_ref <- fit
      pars_ref <- pars
    }
    ref_model <- list(n_obs = n_obs)
    if (!is.null(theta_ref_args$model) && is.list(theta_ref_args$model)) {
      ref_model <- utils::modifyList(ref_model, theta_ref_args$model)
    }
    ref_defaults <- list(
      source = source_ref,
      pars = pars_ref,
      ref = theta_ref,
      mix = mix,
      loglik_fn = loglik_fn,
      model = ref_model
    )
    theta_ref_args_use <- theta_ref_args
    if (!is.null(theta_ref_args_use$model)) theta_ref_args_use$model <- NULL
    ref_call <- utils::modifyList(ref_defaults, theta_ref_args_use)
    theta_hat <- do.call(geoess_theta_ref, ref_call)
  }

  Ipi <- geoess_Ipi_from_mix(mix, theta_hat)

  if (I1_method == "hessian") {
    if (is.null(loglik_fn)) {
      loglik_fn <- .geoess_try_loglik_from_rstan(fit, pars)
    }
    if (is.null(loglik_fn)) {
      stop("loglik_fn is required for Hessian-based I1")
    }
    I1 <- geoess_I1_hessian(theta_hat, loglik_fn, n_obs, ...)
  } else {
    if (is.null(score_fn)) stop("score_fn is required for score-based I1")
    I1 <- geoess_I1_score(theta_hat, score_fn, n_obs, ...)
  }

  ess_out <- if (inherits(mix, "geoess_mix") && method %in% c("trace", "dir")) {
    geoess_ess(I1 = I1, Ipi = mix, method = method, v = v)
  } else {
    geoess_ess(I1 = I1, Ipi = Ipi, method = method, v = v)
  }
  mix_ess <- if (is.list(ess_out)) ess_out else NULL
  if (is.list(ess_out) && !is.null(ess_out$ess)) {
    ess_out <- ess_out$ess
  }

  out <- list(
    ess = ess_out,
    ess_components = if (!is.null(mix_ess)) mix_ess$ess_components else NULL,
    ess_weighted = if (!is.null(mix_ess)) mix_ess$ess_weighted else NULL,
    I1 = I1,
    Ipi = Ipi,
    theta_ref = theta_hat,
    method = method,
    mix = mix,
    metadata = list(
      fit_type = .geoess_fit_type(fit),
      n_obs = n_obs,
      I1_method = I1_method,
      mixture_K = mix$K_selected,
      mixture_criterion = mix$criterion,
      mixture_criterion_value = mix[[mix$criterion]]
    )
  )
  class(out) <- "geoess_result"
  out
}

.geoess_fit_type <- function(fit) {
  if (inherits(fit, "stanfit")) return("rstan")
  if (inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanVB")) return("cmdstanr")
  if (is.matrix(fit)) return("matrix")
  if (is.array(fit)) return("array")
  "unknown"
}

.geoess_try_loglik_from_rstan <- function(fit, pars) {
  if (!inherits(fit, "stanfit")) return(NULL)
  if (!requireNamespace("rstan", quietly = TRUE)) return(NULL)
  if (is.null(fit@stanmodel)) return(NULL)

  all_pars <- names(fit@par_dims)
  if (is.null(all_pars)) return(NULL)
  if (!setequal(pars, all_pars)) {
    warning("rstan log_prob requires full parameter vector; provide loglik_fn")
    return(NULL)
  }

  function(theta) {
    theta_list <- as.list(theta)
    names(theta_list) <- pars
    upars <- rstan::unconstrain_pars(fit@stanmodel, theta_list)
    rstan::log_prob(fit@stanmodel, upars, adjust_transform = TRUE)
  }
}
