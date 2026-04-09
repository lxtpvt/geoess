#' Unit Fisher Information for Theta Block
#'
#' Returns per-observation unit information for the group-level parameters.
#'
#' @param I1_obs Fisher information per observation. Matrix, scalar, list, or function.
#' @param group_index Optional group index for list or function inputs.
#' @param ... Passed to \code{I1_obs} if it is a function.
#'
#' @return Numeric matrix for the theta block.
#' @export
geoess_unit_info_theta <- function(I1_obs, group_index = NULL, ...) {
  if (is.function(I1_obs)) {
    return(.geoess_as_matrix(I1_obs(group_index, ...), "Iunit_theta"))
  }
  if (is.list(I1_obs)) {
    if (is.null(group_index)) stop("group_index is required for list inputs")
    return(.geoess_as_matrix(I1_obs[[group_index]], "Iunit_theta"))
  }
  .geoess_as_matrix(I1_obs, "Iunit_theta")
}

#' Unit Fisher Information for Phi Block
#'
#' Returns per-group unit information for hyperparameters.
#'
#' @param Iphi_group Fisher information per group. Matrix, scalar, list, or function.
#' @param ... Passed to \code{Iphi_group} if it is a function.
#'
#' @return Numeric matrix for the phi block.
#' @export
geoess_unit_info_phi <- function(Iphi_group, ...) {
  if (is.function(Iphi_group)) {
    return(.geoess_as_matrix(Iphi_group(...), "Iunit_phi"))
  }
  .geoess_as_matrix(Iphi_group, "Iunit_phi")
}

#' Block Trace ESS
#'
#' @param Iunit Unit Fisher information for a block.
#' @param Ipi Prior curvature for a block.
#'
#' @return Numeric scalar ESS.
#' @export
geoess_block_trace <- function(Iunit, Ipi) {
  geoess_trace(Iunit, Ipi)
}

#' Block ESS Spectrum
#'
#' @param Iunit Unit Fisher information for a block.
#' @param Ipi Prior curvature for a block.
#'
#' @return A list with eigenvalues and eigenvectors.
#' @export
geoess_block_spectrum <- function(Iunit, Ipi) {
  geoess_spectrum(Iunit, Ipi)
}

.geoess_call_theta_phi <- function(fn, theta, phi, ...) {
  if (!is.function(fn)) return(fn)
  res <- tryCatch(fn(theta, phi, ...), error = function(e) NULL)
  if (!is.null(res)) return(res)
  res <- tryCatch(fn(theta, ...), error = function(e) NULL)
  if (!is.null(res)) return(res)
  fn(phi, ...)
}

.geoess_call_phi_theta <- function(fn, phi, theta, ...) {
  if (!is.function(fn)) return(fn)
  res <- tryCatch(fn(phi, theta, ...), error = function(e) NULL)
  if (!is.null(res)) return(res)
  res <- tryCatch(fn(phi, ...), error = function(e) NULL)
  if (!is.null(res)) return(res)
  fn(theta, ...)
}

#' Level-Specific Local ESS
#'
#' Computes local ESS for a parameter block at a reference point.
#'
#' @param theta_ref Numeric vector for the group-level block.
#' @param phi_ref Numeric vector for the hyperparameter block.
#' @param Iunit_fn Unit information function or matrix.
#' @param Ipi_fn Prior curvature function or matrix.
#' @param block Character. "theta" or "phi".
#' @param method Character. "trace", "eig", or "dir".
#' @param v Direction vector for \code{method = "dir"}.
#' @param conditional Optional list with \code{enable}, \code{idx_a}, \code{idx_b}.
#' @param psd Character. PSD handling for prior curvature.
#' @param symmetrize Logical; if TRUE, symmetrize matrices.
#' @param solver Character. "chol" or "solve".
#' @param ... Passed to user-supplied functions.
#'
#' @return A list with block matrices and ESS summaries.
#' @export
geoess_level_local <- function(theta_ref,
                               phi_ref,
                               Iunit_fn,
                               Ipi_fn,
                               block = c("theta", "phi"),
                               method = c("trace", "eig", "dir"),
                               v = NULL,
                               conditional = list(enable = FALSE, idx_a = NULL, idx_b = NULL),
                               psd = c("none", "clip", "nearPD"),
                               symmetrize = TRUE,
                               solver = c("chol", "solve"),
                               ...) {
  block <- match.arg(block)
  method <- match.arg(method)
  psd <- match.arg(psd)
  solver <- match.arg(solver)

  if (block == "theta") {
    Iunit <- .geoess_call_theta_phi(Iunit_fn, theta_ref, phi_ref, ...)
    Ipi <- .geoess_call_theta_phi(Ipi_fn, theta_ref, phi_ref, ...)
  } else {
    Iunit <- .geoess_call_phi_theta(Iunit_fn, phi_ref, theta_ref, ...)
    Ipi <- .geoess_call_phi_theta(Ipi_fn, phi_ref, theta_ref, ...)
  }

  Iunit <- .geoess_as_matrix(Iunit, "Iunit")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  if (psd != "none") Ipi <- .geoess_psd_adjust(Ipi, psd)

  if (isTRUE(conditional$enable)) {
    res <- geoess_conditional(Iunit, Ipi,
                              idx_a = conditional$idx_a,
                              idx_b = conditional$idx_b,
                              method = if (method == "dir") "trace" else "trace",
                              symmetrize = symmetrize,
                              solver = solver)
    ess <- res$ess_cond
    Iunit <- res$I1_cond
    Ipi <- res$Ipi_cond
  } else if (method == "dir") {
    if (is.null(v)) stop("v is required for method = 'dir'")
    ess <- geoess_dir(Iunit, Ipi, v)
  } else if (method == "eig") {
    ess <- geoess_eig(Iunit, Ipi, symmetrize = symmetrize)
  } else {
    ess <- geoess_trace(Iunit, Ipi, symmetrize = symmetrize)
  }

  list(Iunit = Iunit, Ipi = Ipi, ess = ess, method = method, block = block)
}

#' Level-Specific Global ESS
#'
#' Computes global ESS for a parameter block under a weighting distribution.
#'
#' @param Q Weighting distribution (draws matrix or sampler).
#' @param idx_theta Integer indices for theta block in \code{Q} draws.
#' @param idx_phi Integer indices for phi block in \code{Q} draws.
#' @param block Character. "theta" or "phi".
#' @param Iunit_fn Unit information function or matrix.
#' @param Ipi_fn Prior curvature function or matrix or \code{geoess_mix}.
#' @param method Character. "trace", "eig", or "dir".
#' @param v Direction vector for \code{method = "dir"}.
#' @param conditional Optional list with \code{enable}, \code{idx_a}, \code{idx_b}.
#' @param nsim Integer. Number of draws for sampling.
#' @param seed Integer. Seed for reproducibility.
#' @param weights Optional weights aligned with draws.
#' @param psd Character. PSD handling for prior curvature.
#' @param symmetrize Logical; if TRUE, symmetrize matrices.
#' @param ... Passed to user-supplied functions.
#'
#' @return A \code{geoess_global_result} object.
#' @export
geoess_level_global <- function(Q,
                                idx_theta,
                                idx_phi,
                                block = c("theta", "phi"),
                                Iunit_fn,
                                Ipi_fn,
                                method = c("trace", "eig", "dir"),
                                v = NULL,
                                conditional = list(enable = FALSE, idx_a = NULL, idx_b = NULL),
                                nsim = 2000,
                                seed = 1,
                                weights = NULL,
                                psd = c("none", "clip", "nearPD"),
                                symmetrize = TRUE,
                                ...) {
  block <- match.arg(block)
  method <- match.arg(method)
  psd <- match.arg(psd)

  q <- .geoess_Q_draws(Q, nsim = nsim, seed = seed)
  draws <- q$draws
  n <- nrow(draws)
  d <- ncol(draws)
  if (n == 0) stop("No draws available for Q")

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

  ess_vec <- rep(NA_real_, n)
  eigs <- NULL
  if (method == "eig") {
    block_dim <- if (block == "theta") length(idx_theta) else length(idx_phi)
    eigs <- matrix(NA_real_, nrow = n, ncol = block_dim)
  }

  for (i in seq_len(n)) {
    theta <- draws[i, idx_theta, drop = FALSE]
    phi <- draws[i, idx_phi, drop = FALSE]
    theta <- as.numeric(theta)
    phi <- as.numeric(phi)

    if (block == "theta") {
      Iunit <- .geoess_call_theta_phi(Iunit_fn, theta, phi, ...)
    } else {
      Iunit <- .geoess_call_phi_theta(Iunit_fn, phi, theta, ...)
    }
    Iunit <- .geoess_as_matrix(Iunit, "Iunit")

    if (inherits(Ipi_fn, "geoess_mix")) {
      Ipi_full <- geoess_Ipi_from_mix(Ipi_fn, draws[i, ], psd = psd)
      Ipi <- if (block == "theta") {
        Ipi_full[idx_theta, idx_theta, drop = FALSE]
      } else {
        Ipi_full[idx_phi, idx_phi, drop = FALSE]
      }
    } else {
      if (block == "theta") {
        Ipi <- .geoess_call_theta_phi(Ipi_fn, theta, phi, ...)
      } else {
        Ipi <- .geoess_call_phi_theta(Ipi_fn, phi, theta, ...)
      }
    }
    Ipi <- .geoess_as_matrix(Ipi, "Ipi")
    if (psd != "none") Ipi <- .geoess_psd_adjust(Ipi, psd)

    if (isTRUE(conditional$enable)) {
      res <- geoess_conditional(Iunit, Ipi,
                                idx_a = conditional$idx_a,
                                idx_b = conditional$idx_b,
                                symmetrize = symmetrize)
      ess_vec[i] <- as.numeric(res$ess_cond)
    } else if (method == "dir") {
      ess_vec[i] <- geoess_dir(Iunit, Ipi, v)
    } else if (method == "eig") {
      eigs[i, ] <- geoess_eig(Iunit, Ipi, symmetrize = symmetrize)
    } else {
      ess_vec[i] <- geoess_trace(Iunit, Ipi, symmetrize = symmetrize)
    }
  }

  if (method == "eig") {
    valid <- apply(is.finite(eigs), 1, all)
    if (!all(valid)) {
      warning("Dropping draws with non-finite eigen ESS values")
      eigs <- eigs[valid, , drop = FALSE]
      w_norm <- w_norm[valid]
      w_norm <- w_norm / sum(w_norm)
      n_eff <- .geoess_is_ess(w_norm)
    }
    if (nrow(eigs) == 0) stop("No valid eigen ESS samples")
    estimate <- apply(eigs, 2, .geoess_weighted_mean, w = w_norm)
    mcse <- apply(eigs, 2, function(x) sqrt(.geoess_weighted_var(x, w_norm)) / sqrt(n_eff))
    res <- list(
      type = if (!is.null(q$type)) q$type else "custom",
      method = paste0("level_", block),
      estimate = estimate,
      mcse = mcse,
      n = nrow(eigs),
      weights_summary = if (!is.null(w_all)) .geoess_weight_summary(w_norm) else NULL,
      call = match.call()
    )
    class(res) <- "geoess_global_result"
    return(res)
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

  res <- list(
    type = if (!is.null(q$type)) q$type else "custom",
    method = paste0("level_", block),
    estimate = estimate,
    mcse = mcse,
    n = length(ess_vec),
    weights_summary = if (!is.null(w_all)) .geoess_weight_summary(w_norm) else NULL,
    quantiles = quant,
    call = match.call()
  )
  class(res) <- "geoess_global_result"
  res
}

#' Hierarchical ESS: Level-Specific (Theta and Phi)
#'
#' Computes level-specific ESS for theta and phi blocks, using per-observation
#' and per-group unit information, respectively.
#'
#' @param fit Optional Stan fit.
#' @param draws Optional draws matrix for (theta, phi).
#' @param pars_theta Character vector of theta parameter names.
#' @param pars_phi Character vector of phi parameter names.
#' @param ref Character. "prior_mean", "posterior_mean", "posterior_median",
#'   "posterior_mode", or "custom".
#' @param theta_ref,phi_ref Optional reference points when \code{ref = "custom"}.
#' @param Q Character. "prior", "posterior", or "custom" for global ESS.
#' @param Q_draws Optional draws matrix used when \code{Q = "custom"} or \code{ref = "prior_*"}.
#' @param method Character. "trace", "eig", or "dir".
#' @param conditional Optional list with \code{enable}, \code{idx_a}, \code{idx_b}.
#' @param Iunit_theta_fn,Iunit_phi_fn Unit information functions or matrices.
#' @param Ipi_theta_fn,Ipi_phi_fn Prior curvature functions or matrices.
#' @param Ipi_joint_fn Optional joint curvature function on (theta,phi).
#' @param mix_joint Optional joint mixture for curvature.
#' @param nsim Integer. Number of draws for global ESS.
#' @param seed Integer seed.
#' @param psd Character. PSD handling for Ipi.
#' @param ... Passed to user-supplied functions.
#'
#' @return A list with local and global ESS for theta and phi.
#' @export
geoess_hier_level <- function(fit = NULL,
                              draws = NULL,
                              pars_theta = NULL,
                              pars_phi = NULL,
                              ref = c("prior_mean", "posterior_mean", "posterior_median",
                                      "posterior_mode", "custom"),
                              theta_ref = NULL,
                              phi_ref = NULL,
                              Q = c("prior", "posterior", "custom"),
                              Q_draws = NULL,
                              method = c("trace", "eig", "dir"),
                              conditional = list(enable = FALSE, idx_a = NULL, idx_b = NULL),
                              Iunit_theta_fn = NULL,
                              Iunit_phi_fn = NULL,
                              Ipi_theta_fn = NULL,
                              Ipi_phi_fn = NULL,
                              Ipi_joint_fn = NULL,
                              mix_joint = NULL,
                              nsim = 2000,
                              seed = 1,
                              psd = c("none", "clip", "nearPD"),
                              ...) {
  ref <- match.arg(ref)
  Q <- match.arg(Q)
  method <- match.arg(method)
  psd <- match.arg(psd)

  if (ref == "custom") {
    if (is.null(theta_ref) || is.null(phi_ref)) {
      stop("theta_ref and phi_ref are required for ref = 'custom'")
    }
  } else {
    if (is.null(pars_theta) || is.null(pars_phi)) {
      stop("pars_theta and pars_phi are required for non-custom ref")
    }
    if (is.null(fit) && is.null(draws)) {
      stop("fit or draws are required for non-custom ref")
    }
    if (ref %in% c("posterior_mean", "posterior_median", "posterior_mode")) {
      ref_vals <- geoess_ref_hier(if (!is.null(fit)) fit else draws,
                                  pars_theta = pars_theta, pars_phi = pars_phi,
                                  ref = ref)
      theta_ref <- ref_vals$theta_ref
      phi_ref <- ref_vals$phi_ref
    } else {
      ref_draws <- if (!is.null(Q_draws)) Q_draws else if (!is.null(draws)) draws else NULL
      if (is.null(ref_draws)) stop("Provide Q_draws for prior reference points")
      ref_draws <- as.matrix(ref_draws)
      if (is.null(pars_theta) || is.null(pars_phi)) {
        stop("pars_theta and pars_phi are required for prior reference points")
      }
      if (!is.null(colnames(ref_draws))) {
        ref_draws <- ref_draws[, c(pars_theta, pars_phi), drop = FALSE]
      }
      if (ref == "prior_mean") {
        ref_vec <- colMeans(ref_draws)
      } else {
        ref_vec <- apply(ref_draws, 2, stats::median)
      }
      theta_ref <- ref_vec[seq_len(length(pars_theta))]
      phi_ref <- ref_vec[length(pars_theta) + seq_len(length(pars_phi))]
    }
  }

  if (is.null(Ipi_theta_fn) || is.null(Ipi_phi_fn)) {
    if (!is.null(mix_joint)) {
      Ipi_joint_fn <- function(Theta) geoess_Ipi_from_mix(mix_joint, Theta, psd = psd)
    }
    if (is.function(Ipi_joint_fn)) {
      idx_theta <- seq_len(length(theta_ref))
      idx_phi <- length(theta_ref) + seq_len(length(phi_ref))
      Ipi_theta_fn <- function(theta, phi, ...) {
        Ipi_joint_fn(c(theta, phi))[idx_theta, idx_theta, drop = FALSE]
      }
      Ipi_phi_fn <- function(phi, theta, ...) {
        Ipi_joint_fn(c(theta, phi))[idx_phi, idx_phi, drop = FALSE]
      }
    }
  }

  if (is.null(Iunit_theta_fn) || is.null(Iunit_phi_fn)) {
    stop("Iunit_theta_fn and Iunit_phi_fn are required")
  }
  if (is.null(Ipi_theta_fn) || is.null(Ipi_phi_fn)) {
    stop("Ipi_theta_fn and Ipi_phi_fn are required")
  }

  theta_local <- geoess_level_local(
    theta_ref = theta_ref,
    phi_ref = phi_ref,
    Iunit_fn = Iunit_theta_fn,
    Ipi_fn = Ipi_theta_fn,
    block = "theta",
    method = method,
    conditional = conditional,
    psd = psd,
    ...
  )

  phi_local <- geoess_level_local(
    theta_ref = theta_ref,
    phi_ref = phi_ref,
    Iunit_fn = Iunit_phi_fn,
    Ipi_fn = Ipi_phi_fn,
    block = "phi",
    method = method,
    conditional = conditional,
    psd = psd,
    ...
  )

  global_theta <- NULL
  global_phi <- NULL
  if (Q %in% c("prior", "posterior", "custom")) {
    Q_source <- NULL
    if (Q == "custom") {
      Q_source <- Q_draws
    } else if (Q == "posterior") {
      Q_source <- if (!is.null(fit)) {
        geoess_draws(fit, pars = c(pars_theta, pars_phi), type = "posterior", as = "matrix")
      } else {
        draws
      }
    } else if (Q == "prior") {
      Q_source <- Q_draws
    }
    if (!is.null(Q_source)) {
      Q_source <- as.matrix(Q_source)
      if (!is.null(colnames(Q_source)) && !is.null(pars_theta) && !is.null(pars_phi)) {
        idx_theta <- match(pars_theta, colnames(Q_source))
        idx_phi <- match(pars_phi, colnames(Q_source))
      } else {
        idx_theta <- seq_len(length(theta_ref))
        idx_phi <- length(theta_ref) + seq_len(length(phi_ref))
      }
      global_theta <- geoess_level_global(
        Q = Q_source,
        idx_theta = idx_theta,
        idx_phi = idx_phi,
        block = "theta",
        Iunit_fn = Iunit_theta_fn,
        Ipi_fn = Ipi_theta_fn,
        method = method,
        conditional = conditional,
        nsim = nsim,
        seed = seed,
        psd = psd,
        ...
      )
      global_phi <- geoess_level_global(
        Q = Q_source,
        idx_theta = idx_theta,
        idx_phi = idx_phi,
        block = "phi",
        Iunit_fn = Iunit_phi_fn,
        Ipi_fn = Ipi_phi_fn,
        method = method,
        conditional = conditional,
        nsim = nsim,
        seed = seed,
        psd = psd,
        ...
      )
    }
  }

  list(
    ess_theta = list(local = theta_local, global = global_theta),
    ess_phi = list(local = phi_local, global = global_phi),
    meta = list(ref = ref, Q = Q, method = method, psd = psd)
  )
}

#' Hierarchical Level-Specific ESS from Stan Fits
#'
#' @param fit Stan fit object.
#' @param pars_theta Character vector of theta parameters.
#' @param pars_phi Character vector of phi parameters.
#' @param group_structure Optional list mapping group data.
#' @param loglik_theta_i_fn Optional log-likelihood for one observation.
#' @param loglik_phi_group_fn Optional log-likelihood for one group.
#' @param prior_draws Logical; if TRUE use prior draws when available.
#' @param posterior_draws Logical; if TRUE use posterior draws.
#' @param mixfit Logical; if TRUE fit a joint mixture for Ipi.
#' @param unit_phi Character. "group" (default) or "observation".
#' @param ... Passed to \code{geoess_hier_level}.
#'
#' @return A list with level-specific ESS outputs.
#' @export
geoess_hier_level_from_stan <- function(fit,
                                        pars_theta,
                                        pars_phi,
                                        group_structure = NULL,
                                        loglik_theta_i_fn = NULL,
                                        loglik_phi_group_fn = NULL,
                                        prior_draws = TRUE,
                                        posterior_draws = TRUE,
                                        mixfit = TRUE,
                                        unit_phi = c("group", "observation"),
                                        ...) {
  unit_phi <- match.arg(unit_phi)

  draws_post <- if (posterior_draws) {
    geoess_draws(fit, pars = c(pars_theta, pars_phi), type = "posterior", as = "matrix")
  } else {
    NULL
  }

  mix_joint <- NULL
  if (mixfit && !is.null(draws_post)) {
    mix_joint <- geoess_automixfit(draws_post, family = "mvnorm", K_grid = 1:3, criterion = "BIC")
  }

  Iunit_theta_fn <- NULL
  if (is.function(loglik_theta_i_fn) && !is.null(group_structure)) {
    Iunit_theta_fn <- function(theta, phi, group_index = 1L, ...) {
      y_i <- group_structure[[group_index]]
      if (length(y_i) < 1L) stop("group_structure contains empty group")
      f <- function(th) loglik_theta_i_fn(th, y_i[1], ...)
      -numDeriv::hessian(f, theta)
    }
  }

  Iunit_phi_fn <- NULL
  if (is.function(loglik_phi_group_fn) && !is.null(group_structure)) {
    Iunit_phi_fn <- function(phi, theta, group_index = 1L, ...) {
      y_i <- group_structure[[group_index]]
      f <- function(ph) loglik_phi_group_fn(ph, y_i, ...)
      -numDeriv::hessian(f, phi)
    }
  }

  if (is.null(Iunit_theta_fn) || is.null(Iunit_phi_fn)) {
    stop("Provide loglik_theta_i_fn/loglik_phi_group_fn or unit information functions via ...")
  }

  geoess_hier_level(
    fit = fit,
    draws = draws_post,
    pars_theta = pars_theta,
    pars_phi = pars_phi,
    Iunit_theta_fn = Iunit_theta_fn,
    Iunit_phi_fn = Iunit_phi_fn,
    mix_joint = mix_joint,
    ...
  )
}

#' Hierarchical ESS Dispatcher
#'
#' @param method Character. "level" or "between_borrowed".
#' @param ... Passed to the selected method.
#'
#' @return Hierarchical ESS result.
#' @export
geoess_hier <- function(method = c("level", "between_borrowed"), ...) {
  method <- match.arg(method)
  if (method == "level") {
    geoess_hier_level(...)
  } else {
    geoess_ess_hier(...)
  }
}
