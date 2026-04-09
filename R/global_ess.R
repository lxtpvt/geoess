#' Global ESS as an Expectation of Local ESS
#'
#' Computes global ESS as an expectation of local ESS values under a weighting
#' distribution \eqn{Q(\theta)}.
#'
#' @param Q Distribution for \eqn{\theta}. Can be a draws matrix, list with
#'   \code{draws} and optional \code{weights}, a function \code{rQ(n)},
#'   or a \code{geoess_mix} object.
#' @param I1 Fisher information per observation. Matrix, function, or list
#'   with \code{type = "hessian"}.
#' @param Ipi Prior curvature. Matrix, function, \code{geoess_mix}, or list
#'   with \code{type = "hessian"}.
#' @param method Character. "trace", "dir", or "eig".
#' @param v Direction for \code{method = "dir"}.
#' @param ess_fun Optional custom scalarization \code{ess_fun(I1, Ipi, theta)}.
#' @param psd Character. "none", "clip", or "nearPD" to enforce PSD on Ipi.
#' @param nsim Integer. Number of draws when \code{Q} is a generator or mixture.
#' @param seed Integer seed for reproducibility when sampling.
#' @param weights Optional weights aligned with draws.
#' @param mixture_attribution Character. When \code{Ipi} is a mixture prior,
#'   controls how pointwise local ESS is formed before averaging. The default
#'   \code{"hard_resp"} uses the weighted hard-assignment rule:
#'   \eqn{w_{k^*}\,\mathrm{ESS}_k(\theta)} for a unique max-responsibility
#'   component, with ties averaged across the tied components. \code{"mix"}
#'   uses full-mixture curvature \eqn{-\nabla^2 \log \sum_k w_k \pi_k(\theta)}.
#'   The options \code{"weight"}, \code{"responsibility"}, and \code{"none"}
#'   use \eqn{\sum_k s_k(\theta) I_{\pi_k}(\theta)} with
#'   \eqn{s_k(\theta) = w_k}, \eqn{r_k(\theta)}, and \eqn{1}, respectively.
#' @param return_samples Logical; if TRUE return local ESS samples.
#' @param eig_summary Character. "mean" or "quantiles" for eigen summaries.
#' @param eig_probs Numeric vector of quantile probabilities.
#' @param ... Passed to log-likelihood or log-prior functions.
#'
#' @return A \code{geoess_global_result} object.
#' @export
geoess_global <- function(Q,
                          I1,
                          Ipi,
                          method = c("trace", "dir", "eig"),
                          v = NULL,
                          ess_fun = NULL,
                          psd = c("none", "clip", "nearPD"),
                          nsim = 5000,
                          seed = 1,
                          weights = NULL,
                          mixture_attribution = c("hard_resp", "mix", "weight", "responsibility", "none"),
                          return_samples = FALSE,
                          eig_summary = c("mean", "quantiles"),
                          eig_probs = c(0.05, 0.5, 0.95),
                          ...) {
  method <- match.arg(method)
  psd <- match.arg(psd)
  eig_summary <- match.arg(eig_summary)
  mixture_attribution <- match.arg(mixture_attribution)

  if (!is.null(ess_fun) && method == "eig") {
    stop("ess_fun cannot be used with method = 'eig'")
  }
  if (method == "dir" && is.null(v) && is.null(ess_fun)) {
    stop("v is required for method = 'dir'")
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim <= 0) {
    stop("nsim must be a positive integer")
  }

  q <- .geoess_Q_draws(Q, nsim = nsim, seed = seed)
  draws <- q$draws
  n <- nrow(draws)
  d <- ncol(draws)
  if (n == 0) stop("No draws available for Q")
  if (is.null(d) || d < 1L) stop("Draws must have at least one column")

  w_all <- if (!is.null(q$weights)) q$weights else NULL
  if (!is.null(weights)) {
    if (is.null(w_all)) {
      w_all <- weights
    } else {
      if (length(weights) != length(w_all)) stop("weights length mismatch")
      w_all <- w_all * weights
    }
  }
  w_norm <- .geoess_normalize_weights(w_all, n)
  n_eff <- .geoess_is_ess(w_norm)

  I1_get <- .geoess_I1_getter(I1, n, d, ...)
  Ipi_get <- .geoess_Ipi_getter(Ipi, n, d, psd, mixture_attribution = mixture_attribution, ...)
  I1_scalar <- is.numeric(I1) && length(I1) == 1L

  nonpsd_count <- 0
  ess_samples <- NULL

  if (method == "eig") {
    eigs <- matrix(NA_real_, nrow = n, ncol = d)
    for (i in seq_len(n)) {
      theta <- draws[i, ]
      I1_i <- I1_get(i, theta)
      if (I1_scalar && is.matrix(I1_i) && (nrow(I1_i) == 0 || ncol(I1_i) == 0)) {
        d_loc <- length(theta)
        if (d_loc < 1L) stop("I1 dimension is invalid")
        I1_i <- .geoess_diag_scalar(as.numeric(I1), d_loc)
      }
      Ipi_i <- Ipi_get(i, theta)
      if (any(!is.finite(as.numeric(I1_i))) || any(!is.finite(as.numeric(Ipi_i)))) {
        next
      }
      if (psd == "none" && !.geoess_check_pd(Ipi_i)) nonpsd_count <- nonpsd_count + 1
      eigs[i, ] <- geoess_eigen(I1_i, Ipi_i)
    }

    valid <- apply(is.finite(eigs), 1, all)
    if (!all(valid)) {
      warning("Dropping draws with non-finite eigen ESS")
      eigs <- eigs[valid, , drop = FALSE]
      w_norm <- w_norm[valid]
      w_norm <- w_norm / sum(w_norm)
      n_eff <- .geoess_is_ess(w_norm)
    }
    if (nrow(eigs) == 0) stop("No valid eigen ESS samples")

    estimate <- apply(eigs, 2, .geoess_weighted_mean, w = w_norm)
    mcse <- apply(eigs, 2, function(x) sqrt(.geoess_weighted_var(x, w_norm)) / sqrt(n_eff))
    eig_quant <- NULL
    if (eig_summary == "quantiles") {
      eig_quant <- t(apply(eigs, 2, .geoess_weighted_quantile, w = w_norm, probs = eig_probs))
      rownames(eig_quant) <- paste0("eig", seq_len(ncol(eigs)))
      colnames(eig_quant) <- paste0("q", eig_probs)
    }

    if (return_samples) ess_samples <- eigs

    if (psd == "none" && nonpsd_count > 0) {
      warning("Ipi not PSD for ", nonpsd_count, " draws")
    }

    res <- list(
      type = if (!is.null(q$type)) q$type else "custom",
      method = method,
      estimate = estimate,
      mcse = mcse,
      n = nrow(eigs),
      weights_summary = if (!is.null(w_all)) .geoess_weight_summary(w_norm) else NULL,
      eig_summary = list(summary = eig_summary, probs = eig_probs, quantiles = eig_quant),
      ess_samples = ess_samples,
      meta = list(
        psd = psd,
        nonpsd = nonpsd_count,
        mixture_attribution = if (inherits(Ipi, "geoess_mix") || inherits(Ipi, "betaMix")) mixture_attribution else NULL
      ),
      call = match.call()
    )
    class(res) <- "geoess_global_result"
    return(res)
  }

  ess_vec <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    theta <- draws[i, ]
    I1_i <- I1_get(i, theta)
    if (I1_scalar && is.matrix(I1_i) && (nrow(I1_i) == 0 || ncol(I1_i) == 0)) {
      d_loc <- length(theta)
      if (d_loc < 1L) stop("I1 dimension is invalid")
      I1_i <- .geoess_diag_scalar(as.numeric(I1), d_loc)
    }
    Ipi_i <- Ipi_get(i, theta)
    if (any(!is.finite(as.numeric(I1_i))) || any(!is.finite(as.numeric(Ipi_i)))) {
      next
    }
    if (psd == "none" && !.geoess_check_pd(Ipi_i)) nonpsd_count <- nonpsd_count + 1

    if (!is.null(ess_fun)) {
      ess_vec[i] <- ess_fun(I1_i, Ipi_i, theta)
    } else if (method == "trace") {
      ess_vec[i] <- geoess_trace(I1_i, Ipi_i)
    } else {
      ess_vec[i] <- geoess_dir(I1_i, Ipi_i, v)
    }
  }

  valid <- is.finite(ess_vec)
  if (!all(valid)) {
    warning("Dropping draws with non-finite ESS values")
    ess_vec <- ess_vec[valid]
    w_norm <- w_norm[valid]
    w_norm <- w_norm / sum(w_norm)
    n_eff <- .geoess_is_ess(w_norm)
  }
  if (length(ess_vec) == 0) stop("No valid ESS samples")

  estimate <- .geoess_weighted_mean(ess_vec, w_norm)
  mcse <- sqrt(.geoess_weighted_var(ess_vec, w_norm)) / sqrt(n_eff)
  quant_probs <- c(0.05, 0.5, 0.95)
  quant <- .geoess_weighted_quantile(ess_vec, w_norm, probs = quant_probs)
  names(quant) <- paste0("q", quant_probs)

  if (return_samples) ess_samples <- ess_vec

  if (psd == "none" && nonpsd_count > 0) {
    warning("Ipi not PSD for ", nonpsd_count, " draws")
  }

  res <- list(
    type = if (!is.null(q$type)) q$type else "custom",
    method = method,
    estimate = estimate,
    mcse = mcse,
    n = length(ess_vec),
    weights_summary = if (!is.null(w_all)) .geoess_weight_summary(w_norm) else NULL,
    quantiles = quant,
    ess_samples = ess_samples,
    meta = list(
      psd = psd,
      nonpsd = nonpsd_count,
      mixture_attribution = if (inherits(Ipi, "geoess_mix") || inherits(Ipi, "betaMix")) mixture_attribution else NULL
    ),
    call = match.call()
  )
  class(res) <- "geoess_global_result"
  res
}

#' Global ESS with Prior Weights
#'
#' @param prior Prior representation: draws matrix, \code{geoess_mix}, or Stan fit.
#' @param I1 Fisher information per observation.
#' @param Ipi Prior curvature or NULL to derive from prior draws.
#' @param pars Optional parameter names for Stan fits.
#' @param ipi_source Character. "mixture", "num_hessian", or "fixed".
#' @param mixture_args Optional list passed to \code{geoess_automixfit}.
#' @param ... Passed to \code{geoess_global}.
#'
#' @return A \code{geoess_global_result}.
#' @export
geoess_global_prior <- function(prior,
                                I1,
                                Ipi = NULL,
                                pars = NULL,
                                ipi_source = c("mixture", "num_hessian", "fixed"),
                                mixture_args = list(),
                                ...) {
  ipi_source <- match.arg(ipi_source)
  q <- .geoess_prior_draws(prior, pars = pars)

  dots <- list(...)
  logprior_fn <- dots$logprior_fn
  if (!is.null(dots$Ipi)) dots$Ipi <- NULL
  if (!is.null(dots$logprior_fn)) dots$logprior_fn <- NULL

  if (is.null(Ipi)) {
    if (ipi_source == "mixture") {
      mix <- if (inherits(prior, "geoess_mix")) prior else do.call(geoess_automixfit, c(list(q), mixture_args))
      Ipi <- mix
    } else if (ipi_source == "num_hessian") {
      if (is.null(logprior_fn)) stop("logprior_fn is required for ipi_source = 'num_hessian'")
      Ipi <- list(type = "hessian", logprior_fn = logprior_fn)
    } else {
      stop("Ipi must be provided when ipi_source = 'fixed'")
    }
  }

  Q_arg <- if (inherits(q, "geoess_mix")) q else list(draws = q, type = "prior")
  res <- do.call(geoess_global, c(list(Q = Q_arg, I1 = I1, Ipi = Ipi), dots))
  res$type <- "prior"
  res
}

#' Global ESS with Posterior Weights
#'
#' @param posterior Posterior representation: draws matrix, \code{geoess_mix}, or Stan fit.
#' @param I1 Fisher information per observation.
#' @param Ipi Prior curvature or NULL to derive from posterior draws.
#' @param pars Optional parameter names for Stan fits.
#' @param ipi_source Character. "mixture", "num_hessian", or "fixed".
#' @param mixture_args Optional list passed to \code{geoess_automixfit}.
#' @param ... Passed to \code{geoess_global}.
#'
#' @return A \code{geoess_global_result}.
#' @export
geoess_global_posterior <- function(posterior,
                                    I1,
                                    Ipi = NULL,
                                    pars = NULL,
                                    ipi_source = c("mixture", "num_hessian", "fixed"),
                                    mixture_args = list(),
                                    ...) {
  ipi_source <- match.arg(ipi_source)
  q <- .geoess_posterior_draws(posterior, pars = pars)

  dots <- list(...)
  logprior_fn <- dots$logprior_fn
  if (!is.null(dots$Ipi)) dots$Ipi <- NULL
  if (!is.null(dots$logprior_fn)) dots$logprior_fn <- NULL

  if (is.null(Ipi)) {
    if (ipi_source == "mixture") {
      mix <- if (inherits(posterior, "geoess_mix")) posterior else do.call(geoess_automixfit, c(list(q), mixture_args))
      Ipi <- mix
    } else if (ipi_source == "num_hessian") {
      if (is.null(logprior_fn)) stop("logprior_fn is required for ipi_source = 'num_hessian'")
      Ipi <- list(type = "hessian", logprior_fn = logprior_fn)
    } else {
      stop("Ipi must be provided when ipi_source = 'fixed'")
    }
  }

  Q_arg <- if (inherits(q, "geoess_mix")) q else list(draws = q, type = "posterior")
  res <- do.call(geoess_global, c(list(Q = Q_arg, I1 = I1, Ipi = Ipi), dots))
  res$type <- "posterior"
  res
}

#' Predictive-Weighted Global ESS
#'
#' @param prior Prior representation: draws matrix, \code{geoess_mix}, or Stan fit.
#' @param I1 Fisher information per observation.
#' @param yrep Optional fixed predictive data.
#' @param loglik_point Function \code{log p(yrep | theta)}.
#' @param r_yrep Optional generator for predictive data.
#' @param weight_mode Character. "analytic" or "mc".
#' @param pars Optional parameter names for Stan fits.
#' @param ipi_source Character. "mixture", "num_hessian", or "fixed".
#' @param mixture_args Optional list passed to \code{geoess_automixfit}.
#' @param ... Passed to \code{geoess_global}.
#'
#' @return A \code{geoess_global_result}.
#' @export
geoess_global_predictive <- function(prior,
                                     I1,
                                     yrep = NULL,
                                     loglik_point = NULL,
                                     r_yrep = NULL,
                                     weight_mode = c("analytic", "mc"),
                                     pars = NULL,
                                     ipi_source = c("mixture", "num_hessian", "fixed"),
                                     mixture_args = list(),
                                     ...) {
  weight_mode <- match.arg(weight_mode)
  ipi_source <- match.arg(ipi_source)
  if (is.null(loglik_point)) stop("loglik_point is required")

  q <- .geoess_prior_draws(prior, pars = pars)
  q_draws <- if (inherits(q, "geoess_mix")) {
    nsim_use <- list(...)$nsim
    if (is.null(nsim_use)) nsim_use <- 5000
    .geoess_rmix_geoess(q, nsim_use)
  } else {
    q
  }

  if (is.null(yrep)) {
    if (is.null(r_yrep)) stop("Provide yrep or r_yrep")
    yrep <- r_yrep(1)
  }

  dots <- list(...)
  logprior_fn <- dots$logprior_fn
  if (!is.null(dots$logprior_fn)) dots$logprior_fn <- NULL

  weight_args <- .geoess_filter_args(loglik_point, dots)
  if (!is.null(weight_args$logprior_fn)) weight_args$logprior_fn <- NULL
  if (!is.null(weight_args$Ipi)) weight_args$Ipi <- NULL

  logw <- .geoess_loglik_point_eval(loglik_point, q_draws, yrep,
                                    weight_mode = weight_mode, weight_args)
  w <- .geoess_normalize_logw(logw)$weights

  Ipi <- NULL
  if (ipi_source == "mixture") {
    mix <- if (inherits(prior, "geoess_mix")) prior else do.call(geoess_automixfit, c(list(q_draws), mixture_args))
    Ipi <- mix
    if (!is.null(dots$Ipi)) dots$Ipi <- NULL
  } else if (ipi_source == "num_hessian") {
    if (is.null(logprior_fn)) stop("logprior_fn is required for ipi_source = 'num_hessian'")
    Ipi <- list(type = "hessian", logprior_fn = logprior_fn)
    if (!is.null(dots$Ipi)) dots$Ipi <- NULL
  } else {
    Ipi <- dots$Ipi
    if (is.null(Ipi)) stop("Ipi must be provided when ipi_source = 'fixed'")
    dots$Ipi <- NULL
  }

  res <- do.call(geoess_global, c(list(Q = list(draws = q_draws, weights = w, type = "predictive"),
                                       I1 = I1, Ipi = Ipi), dots))
  res$type <- "predictive"
  res
}

#' Design-Based Global ESS
#'
#' @param rQ Function \code{rQ(n)} returning draws from a design distribution.
#' @param I1 Fisher information per observation.
#' @param Ipi Prior curvature.
#' @param ... Passed to \code{geoess_global}.
#'
#' @return A \code{geoess_global_result}.
#' @export
geoess_global_design <- function(rQ, I1, Ipi, ...) {
  if (!is.function(rQ)) stop("rQ must be a function")
  res <- geoess_global(Q = rQ, I1 = I1, Ipi = Ipi, ...)
  res$type <- "design"
  res
}

.geoess_Q_draws <- function(Q, nsim, seed) {
  if (is.matrix(Q) || is.array(Q)) {
    draws <- as.matrix(Q)
    return(list(draws = draws, weights = NULL, type = "custom"))
  }

  if (inherits(Q, "geoess_spec") || (is.list(Q) && !is.null(Q$prior_source))) {
    type <- if (!is.null(Q$type)) Q$type else if (!is.null(Q$posterior_source)) "posterior" else "prior"
    source <- if (type == "prior") Q$prior_source else Q$posterior_source
    if (is.null(source)) source <- Q$source
    if (is.null(source)) stop("geoess_spec does not contain a source")
    draws <- geoess_draws(source, pars = Q$pars, type = type, as = "matrix")
    weights <- Q$weights
    return(list(draws = draws, weights = weights, type = type))
  }

  if (inherits(Q, "geoess_mix")) {
    if (!is.null(seed)) set.seed(seed)
    draws <- .geoess_rmix_geoess(Q, nsim)
    return(list(draws = draws, weights = NULL, type = "custom"))
  }

  if (is.function(Q)) {
    if (!is.null(seed)) set.seed(seed)
    draws <- Q(nsim)
    return(list(draws = as.matrix(draws), weights = NULL, type = "custom"))
  }

  if (is.list(Q) && !is.null(Q$draws)) {
    draws <- as.matrix(Q$draws)
    weights <- Q$weights
    type <- if (!is.null(Q$type)) Q$type else "custom"
    return(list(draws = draws, weights = weights, type = type))
  }

  if (is.list(Q) && !is.null(Q$rQ)) {
    if (!is.null(seed)) set.seed(seed)
    draws <- Q$rQ(nsim)
    weights <- Q$weights
    type <- if (!is.null(Q$type)) Q$type else "custom"
    return(list(draws = as.matrix(draws), weights = weights, type = type))
  }

  stop("Unsupported Q type")
}

.geoess_I1_getter <- function(I1, n, d, ...) {
  if (is.numeric(I1) && length(I1) == 1L) {
    if (!is.finite(I1)) stop("I1 must be finite")
    I1_val <- as.numeric(I1)
    return(function(i, theta) {
      d_loc <- if (!is.null(d) && d > 0L) d else length(theta)
      if (is.null(d_loc) || d_loc < 1L) stop("I1 dimension is invalid")
      .geoess_diag_scalar(I1_val, d_loc)
    })
  }

  if (is.matrix(I1) || inherits(I1, "Matrix")) {
    I1_const <- .geoess_as_matrix(I1, "I1")
    .geoess_check_square(I1_const, "I1")
    return(function(i, theta) I1_const)
  }

  if (is.function(I1)) {
    return(function(i, theta) {
      I1_mat <- .geoess_as_matrix(I1(theta, ...), "I1")
      .geoess_check_square(I1_mat, "I1")
      I1_mat
    })
  }

  if (is.list(I1) && !is.null(I1$type) && I1$type == "hessian") {
    if (is.null(I1$loglik_fn)) stop("I1 loglik_fn is required")
    if (is.null(I1$n_obs)) stop("I1 n_obs is required")
    return(function(i, theta) geoess_I1_hessian(theta, I1$loglik_fn, I1$n_obs, ...))
  }

  if (is.list(I1) && length(I1) == n) {
    return(function(i, theta) {
      I1_mat <- .geoess_as_matrix(I1[[i]], "I1")
      .geoess_check_square(I1_mat, "I1")
      I1_mat
    })
  }

  if (is.array(I1) && length(dim(I1)) == 3L && dim(I1)[3] == n) {
    return(function(i, theta) {
      I1_mat <- .geoess_as_matrix(I1[, , i], "I1")
      .geoess_check_square(I1_mat, "I1")
      I1_mat
    })
  }

  stop("Unsupported I1 specification")
}

.geoess_Ipi_getter <- function(Ipi, n, d, psd, mixture_attribution = "hard_resp", ...) {
  if (inherits(Ipi, "geoess_mix") || inherits(Ipi, "betaMix")) {
    return(function(i, theta) {
      .geoess_mix_Ipi_global(
        mix = Ipi,
        theta = theta,
        psd = psd,
        mixture_attribution = mixture_attribution,
        ...
      )
    })
  }

  if (is.numeric(Ipi) && length(Ipi) == 1L) {
    if (!is.finite(Ipi)) stop("Ipi must be finite")
    Ipi_val <- as.numeric(Ipi)
    return(function(i, theta) {
      d_loc <- if (!is.null(d) && d > 0L) d else length(theta)
      if (is.null(d_loc) || d_loc < 1L) stop("Ipi dimension is invalid")
      Ipi_const <- .geoess_diag_scalar(Ipi_val, d_loc)
      .geoess_psd_adjust(Ipi_const, psd)
    })
  }

  if (is.matrix(Ipi) || inherits(Ipi, "Matrix")) {
    Ipi_const <- .geoess_as_matrix(Ipi, "Ipi")
    .geoess_check_square(Ipi_const, "Ipi")
    Ipi_const <- .geoess_psd_adjust(Ipi_const, psd)
    return(function(i, theta) Ipi_const)
  }

  if (is.function(Ipi)) {
    return(function(i, theta) {
      Ipi_mat <- .geoess_as_matrix(Ipi(theta, ...), "Ipi")
      .geoess_check_square(Ipi_mat, "Ipi")
      .geoess_psd_adjust(Ipi_mat, psd)
    })
  }

  if (is.list(Ipi) && !is.null(Ipi$type) && Ipi$type == "hessian") {
    if (is.null(Ipi$logprior_fn)) stop("Ipi logprior_fn is required")
    return(function(i, theta) {
      H <- -numDeriv::hessian(Ipi$logprior_fn, theta, ...)
      H <- .geoess_symmetrize(H)
      .geoess_psd_adjust(H, psd)
    })
  }

  if (is.list(Ipi) && length(Ipi) == n) {
    return(function(i, theta) {
      Ipi_mat <- .geoess_as_matrix(Ipi[[i]], "Ipi")
      .geoess_check_square(Ipi_mat, "Ipi")
      .geoess_psd_adjust(Ipi_mat, psd)
    })
  }

  if (is.array(Ipi) && length(dim(Ipi)) == 3L && dim(Ipi)[3] == n) {
    return(function(i, theta) {
      Ipi_mat <- .geoess_as_matrix(Ipi[, , i], "Ipi")
      .geoess_check_square(Ipi_mat, "Ipi")
      .geoess_psd_adjust(Ipi_mat, psd)
    })
  }

  stop("Unsupported Ipi specification")
}

.geoess_mix_Ipi_global <- function(mix,
                                   theta,
                                   psd = c("none", "clip", "nearPD"),
                                   mixture_attribution = c("hard_resp", "mix", "weight", "responsibility", "none"),
                                   use_cache = TRUE,
                                   symmetrize = TRUE,
                                   ...) {
  psd <- match.arg(psd)
  mixture_attribution <- match.arg(mixture_attribution)

  if (mixture_attribution == "mix") {
    return(.geoess_psd_adjust(
      .geoess_mix_Ipi_full(
        mix = mix,
        theta = theta,
        use_cache = use_cache,
        symmetrize = symmetrize,
        ...
      ),
      psd
    ))
  }

  w <- .geoess_mix_weights(mix)
  K <- length(w)
  theta <- as.numeric(theta)

  scales <- switch(
    mixture_attribution,
    weight = w,
    responsibility = .geoess_mix_responsibility(mix, theta, use_cache = use_cache, ...),
    hard_resp = .geoess_hard_resp_scales(
      .geoess_mix_responsibility(mix, theta, use_cache = use_cache, ...),
      w
    ),
    none = rep(1, K)
  )

  Ipi_sum <- NULL
  for (k in seq_len(K)) {
    if (!is.finite(scales[k]) || scales[k] == 0) next
    Ipi_k <- .geoess_mix_component_Ipi(
      mix = mix,
      k = k,
      theta = theta,
      use_cache = use_cache,
      symmetrize = symmetrize,
      ...
    )
    if (is.null(Ipi_sum)) {
      Ipi_sum <- scales[k] * Ipi_k
    } else {
      Ipi_sum <- Ipi_sum + scales[k] * Ipi_k
    }
  }

  if (is.null(Ipi_sum)) {
    d_loc <- .geoess_mix_dim(mix)
    Ipi_sum <- matrix(0, nrow = d_loc, ncol = d_loc)
  }

  .geoess_psd_adjust(Ipi_sum, psd)
}

.geoess_rmix_geoess <- function(mix, n) {
  if (!inherits(mix, "geoess_mix")) stop("mix must be a geoess_mix object")
  if (is.null(mix$family) || mix$family != "mvnorm") {
    stop("Only mvnorm mixtures are supported")
  }
  w <- as.numeric(mix$weights)
  w <- w / sum(w)
  K <- length(w)
  d <- length(mix$mu[[1]])
  comp <- sample(seq_len(K), size = n, replace = TRUE, prob = w)
  out <- matrix(NA_real_, nrow = n, ncol = d)
  for (k in seq_len(K)) {
    idx <- which(comp == k)
    if (!length(idx)) next
    mu_k <- as.numeric(mix$mu[[k]])
    Sig_k <- mix$Sigma[[k]]
    out[idx, ] <- mvtnorm::rmvnorm(length(idx), mean = mu_k, sigma = Sig_k)
  }
  out
}

.geoess_prior_draws <- function(prior, pars = NULL) {
  if (inherits(prior, "geoess_mix")) return(prior)
  if (inherits(prior, "stanfit") || inherits(prior, "CmdStanMCMC") || inherits(prior, "CmdStanVB")) {
    return(geoess_draws(prior, pars = pars, type = "prior", as = "matrix"))
  }
  if (is.matrix(prior) || is.array(prior)) return(as.matrix(prior))
  stop("Unsupported prior input")
}

.geoess_posterior_draws <- function(posterior, pars = NULL) {
  if (inherits(posterior, "geoess_mix")) return(posterior)
  if (inherits(posterior, "stanfit") || inherits(posterior, "CmdStanMCMC") || inherits(posterior, "CmdStanVB")) {
    return(geoess_draws(posterior, pars = pars, type = "posterior", as = "matrix"))
  }
  if (is.matrix(posterior) || is.array(posterior)) return(as.matrix(posterior))
  stop("Unsupported posterior input")
}

.geoess_loglik_point_eval <- function(loglik_point, draws, yrep, weight_mode, args = list()) {
  draws <- as.matrix(draws)
  n <- nrow(draws)

  val <- tryCatch(do.call(loglik_point, c(list(draws, yrep), args)),
                  error = function(e) NULL)
  if (!is.null(val)) {
    if (is.numeric(val) && length(val) == n) return(as.numeric(val))
    if (is.matrix(val)) {
      if (nrow(val) == n) return(rowSums(val))
      if (ncol(val) == n) return(colSums(val))
    }
  }

  logw <- numeric(n)
  for (i in seq_len(n)) {
    logw[i] <- do.call(loglik_point, c(list(draws[i, ], yrep), args))
  }
  logw
}

.geoess_filter_args <- function(fun, args) {
  if (!is.function(fun) || length(args) == 0) return(args)
  fmls <- names(formals(fun))
  if (is.null(fmls) || "..." %in% fmls) return(args)
  args[intersect(names(args), fmls)]
}
