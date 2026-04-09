#' Schur Complement for Hierarchical Curvature
#'
#' Computes conditional curvature for a target block given a nuisance block.
#'
#' @param Htt Numeric matrix. Hessian block for target parameter.
#' @param Htp Numeric matrix. Hessian block for interaction between target and nuisance.
#' @param Hpp Numeric matrix. Hessian block for nuisance parameter.
#' @param symmetrize Logical; if TRUE, symmetrize outputs.
#' @param solver Character. "chol" (preferred) or "solve".
#' @param return_components Logical; if TRUE, return a list with effective,
#'   between-group, and borrowed curvature components.
#'
#' @return When \code{return_components = TRUE}, a list with elements
#'   \code{Ipi_eff}, \code{Ipi_between}, and \code{Ipi_borrowed}. Otherwise,
#'   returns \code{Ipi_eff}.
#'
#' @examples
#' Htt <- matrix(2, 1, 1)
#' Htp <- matrix(1, 1, 1)
#' Hpp <- matrix(3, 1, 1)
#' schur_theta_given_phi(Htt, Htp, Hpp)
#' @export
schur_theta_given_phi <- function(Htt, Htp, Hpp,
                                  symmetrize = TRUE,
                                  solver = c("chol", "solve"),
                                  return_components = TRUE) {
  solver <- match.arg(solver)
  Htt <- .geoess_as_matrix(Htt, "Htt")
  Htp <- .geoess_as_matrix(Htp, "Htp")
  Hpp <- .geoess_as_matrix(Hpp, "Hpp")
  .geoess_check_square(Htt, "Htt")
  .geoess_check_square(Hpp, "Hpp")

  if (nrow(Htp) != nrow(Htt) || ncol(Htp) != nrow(Hpp)) {
    stop("Dimension mismatch between Htt, Htp, and Hpp")
  }

  Hpt <- t(Htp)
  tmp <- NULL
  if (solver == "chol") {
    chol_fit <- .geoess_try_chol(as.matrix(Hpp))
    if (!is.null(chol_fit$R)) {
      R <- chol_fit$R
      tmp <- backsolve(R, forwardsolve(t(R), Hpt))
    } else {
      warning("Hpp is not positive definite; using solve fallback")
    }
  }
  if (is.null(tmp)) {
    tmp <- .geoess_solve(Hpp, Hpt)
  }

  Ipi_eff <- Htt - Htp %*% tmp
  Ipi_between <- Htt
  Ipi_borrowed <- -(Htp %*% tmp)

  if (symmetrize) {
    Ipi_eff <- .geoess_symmetrize(Ipi_eff)
    Ipi_between <- .geoess_symmetrize(Ipi_between)
    Ipi_borrowed <- .geoess_symmetrize(Ipi_borrowed)
  }

  if (return_components) {
    return(list(
      Ipi_eff = Ipi_eff,
      Ipi_between = Ipi_between,
      Ipi_borrowed = Ipi_borrowed
    ))
  }

  Ipi_eff
}

#' Hierarchical Curvature via Schur Complement
#'
#' @param Htt Numeric matrix. Curvature block for \eqn{\theta}.
#' @param Htp Numeric matrix. Cross curvature block.
#' @param Hpp Numeric matrix. Curvature block for \eqn{\phi}.
#' @param symmetrize Logical; if TRUE, symmetrize outputs.
#' @param solver Character. "chol" or "solve".
#' @param return_components Logical; if TRUE, return all blocks.
#'
#' @return A list with effective/between/borrowed curvature blocks.
#' @export
geoess_Ipi_hier <- function(Htt, Htp, Hpp,
                            symmetrize = TRUE,
                            solver = c("chol", "solve"),
                            return_components = TRUE) {
  schur_theta_given_phi(Htt, Htp, Hpp,
                        symmetrize = symmetrize,
                        solver = solver,
                        return_components = return_components)
}

.geoess_split_blocks <- function(H, d_theta, d_phi) {
  idx_theta <- seq_len(d_theta)
  idx_phi <- d_theta + seq_len(d_phi)
  list(
    Htt = H[idx_theta, idx_theta, drop = FALSE],
    Htp = H[idx_theta, idx_phi, drop = FALSE],
    Hpp = H[idx_phi, idx_phi, drop = FALSE]
  )
}

#' Hierarchical Curvature from Functions
#'
#' @param theta_ref Numeric vector for \eqn{\theta^*}.
#' @param phi_ref Numeric vector for \eqn{\phi^*}.
#' @param Htt_fn Optional function returning Htt.
#' @param Htp_fn Optional function returning Htp.
#' @param Hpp_fn Optional function returning Hpp.
#' @param logprior_joint Optional joint log-prior function.
#' @param ... Passed to user functions.
#'
#' @return A list with curvature blocks and effective components.
#' @export
geoess_Ipi_hier_from_fn <- function(theta_ref, phi_ref,
                                    Htt_fn = NULL,
                                    Htp_fn = NULL,
                                    Hpp_fn = NULL,
                                    logprior_joint = NULL,
                                    ...) {
  theta_ref <- as.numeric(theta_ref)
  phi_ref <- as.numeric(phi_ref)
  d_theta <- length(theta_ref)
  d_phi <- length(phi_ref)

  if (!is.null(Htt_fn) && !is.null(Htp_fn) && !is.null(Hpp_fn)) {
    Htt <- Htt_fn(theta_ref, phi_ref, ...)
    Htp <- Htp_fn(theta_ref, phi_ref, ...)
    Hpp <- Hpp_fn(theta_ref, phi_ref, ...)
  } else {
    if (!is.function(logprior_joint)) {
      stop("Provide Htt_fn/Htp_fn/Hpp_fn or logprior_joint")
    }
    x <- c(theta_ref, phi_ref)
    H_joint <- -numDeriv::hessian(logprior_joint, x, ...)
    H_joint <- .geoess_symmetrize(H_joint)
    blocks <- .geoess_split_blocks(H_joint, d_theta, d_phi)
    Htt <- blocks$Htt
    Htp <- blocks$Htp
    Hpp <- blocks$Hpp
  }

  res <- geoess_Ipi_hier(Htt, Htp, Hpp, return_components = TRUE)
  res$Htt <- .geoess_as_matrix(Htt, "Htt")
  res$Htp <- .geoess_as_matrix(Htp, "Htp")
  res$Hpp <- .geoess_as_matrix(Hpp, "Hpp")
  res
}

#' Hierarchical Curvature from Joint Mixture
#'
#' @param mix_joint A \code{geoess_mix} object on joint (theta, phi).
#' @param theta_ref Numeric vector for \eqn{\theta^*}.
#' @param phi_ref Numeric vector for \eqn{\phi^*}.
#' @param psd Character. Passed to \code{geoess_Ipi_from_mix}.
#' @param ... Passed to mixture calculus.
#'
#' @return A list with curvature blocks and effective components.
#' @export
geoess_Ipi_hier_from_mix <- function(mix_joint, theta_ref, phi_ref,
                                     psd = c("none", "clip", "nearPD"),
                                     ...) {
  psd <- match.arg(psd)
  theta_ref <- as.numeric(theta_ref)
  phi_ref <- as.numeric(phi_ref)
  d_theta <- length(theta_ref)
  d_phi <- length(phi_ref)
  theta_phi <- c(theta_ref, phi_ref)

  Ipi_joint <- geoess_Ipi_from_mix(mix_joint, theta_phi, psd = psd, ...)
  blocks <- .geoess_split_blocks(Ipi_joint, d_theta, d_phi)

  res <- geoess_Ipi_hier(blocks$Htt, blocks$Htp, blocks$Hpp, return_components = TRUE)
  res$Htt <- blocks$Htt
  res$Htp <- blocks$Htp
  res$Hpp <- blocks$Hpp
  res
}

#' Reference Point for Hierarchical Parameters
#'
#' @param fit_or_draws Stan fit or draws matrix.
#' @param pars_theta Character vector of theta parameter names.
#' @param pars_phi Character vector of phi parameter names.
#' @param ref Character. One of "posterior_mean", "posterior_median", or "posterior_mode".
#' @param ... Passed to Stan optimization if available.
#'
#' @return A list with \code{theta_ref}, \code{phi_ref}, and \code{ref}.
#' @export
geoess_ref_hier <- function(fit_or_draws, pars_theta, pars_phi,
                            ref = c("posterior_mean", "posterior_median", "posterior_mode"),
                            ...) {
  ref <- match.arg(ref)

  if (inherits(fit_or_draws, "stanfit") ||
      inherits(fit_or_draws, "CmdStanMCMC") ||
      inherits(fit_or_draws, "CmdStanVB")) {
    draws <- geoess_draws(fit_or_draws, pars = c(pars_theta, pars_phi),
                          type = "posterior", as = "matrix")
  } else if (is.matrix(fit_or_draws) || is.array(fit_or_draws)) {
    draws <- as.matrix(fit_or_draws)
  } else {
    stop("fit_or_draws must be a Stan fit or draws matrix")
  }

  if (!is.null(colnames(draws))) {
    missing <- setdiff(c(pars_theta, pars_phi), colnames(draws))
    if (length(missing) > 0) stop("Missing parameters in draws: ", paste(missing, collapse = ", "))
    draws <- draws[, c(pars_theta, pars_phi), drop = FALSE]
  }

  d_theta <- length(pars_theta)
  d_phi <- length(pars_phi)
  if (ncol(draws) != d_theta + d_phi) {
    stop("Dimension mismatch between draws and parameter names")
  }

  if (ref == "posterior_mean") {
    ref_vec <- colMeans(draws)
  } else if (ref == "posterior_median") {
    ref_vec <- apply(draws, 2, stats::median)
  } else {
    if (inherits(fit_or_draws, "stanfit") ||
        inherits(fit_or_draws, "CmdStanMCMC") ||
        inherits(fit_or_draws, "CmdStanVB")) {
      map_val <- .geoess_theta_ref_map(fit_or_draws, pars = c(pars_theta, pars_phi), ...)
      if (!is.null(map_val)) {
        ref_vec <- map_val
      } else {
        warning("MAP not available; falling back to posterior mean")
        ref_vec <- colMeans(draws)
      }
    } else {
      warning("MAP not available for draws; falling back to posterior mean")
      ref_vec <- colMeans(draws)
    }
  }

  list(
    theta_ref = as.numeric(ref_vec[seq_len(d_theta)]),
    phi_ref = as.numeric(ref_vec[d_theta + seq_len(d_phi)]),
    ref = ref
  )
}

#' Hierarchical ESS from Curvature Blocks
#'
#' @param I1 Fisher information matrix or function \code{I1(theta)}.
#' @param Ipi_blocks Optional list with \code{Htt}, \code{Htp}, \code{Hpp}.
#' @param Htt,Htp,Hpp Optional curvature blocks if \code{Ipi_blocks} is NULL.
#' @param theta_ref,phi_ref Reference points for evaluation.
#' @param method Character. "trace", "dir", or "eig".
#' @param v Direction for \code{method = "dir"}.
#' @param psd Character. "none", "clip", or "nearPD" applied to Ipi_eff only.
#' @param return_components Logical; if TRUE, return components.
#' @param ... Passed to \code{I1} if it is a function.
#'
#' @return A list with ESS totals and component decompositions.
#' @export
geoess_ess_hier <- function(I1,
                            Ipi_blocks = NULL,
                            Htt = NULL, Htp = NULL, Hpp = NULL,
                            theta_ref,
                            phi_ref = NULL,
                            method = c("trace", "dir", "eig"),
                            v = NULL,
                            psd = c("none", "clip", "nearPD"),
                            return_components = TRUE,
                            ...) {
  method <- match.arg(method)
  psd <- match.arg(psd)

  Ipi_eff <- NULL
  Ipi_between <- NULL
  Ipi_borrowed <- NULL

  if (!is.null(Ipi_blocks)) {
    if (!is.null(Ipi_blocks$Ipi_eff)) {
      Ipi_eff <- .geoess_as_matrix(Ipi_blocks$Ipi_eff, "Ipi_eff")
      Ipi_between <- .geoess_as_matrix(Ipi_blocks$Ipi_between, "Ipi_between")
      Ipi_borrowed <- .geoess_as_matrix(Ipi_blocks$Ipi_borrowed, "Ipi_borrowed")
      Ipi_eff <- .geoess_symmetrize(Ipi_eff)
      Ipi_between <- .geoess_symmetrize(Ipi_between)
      Ipi_borrowed <- .geoess_symmetrize(Ipi_borrowed)
    } else {
      Htt <- Ipi_blocks$Htt
      Htp <- Ipi_blocks$Htp
      Hpp <- Ipi_blocks$Hpp
      if (is.null(Htt) || is.null(Htp) || is.null(Hpp)) {
        stop("Ipi_blocks must contain Htt/Htp/Hpp or Ipi_eff/Ipi_between/Ipi_borrowed")
      }
      blocks <- geoess_Ipi_hier(Htt, Htp, Hpp, return_components = TRUE)
      Ipi_eff <- blocks$Ipi_eff
      Ipi_between <- blocks$Ipi_between
      Ipi_borrowed <- blocks$Ipi_borrowed
    }
  } else {
    if (is.null(Htt) || is.null(Htp) || is.null(Hpp)) {
      stop("Provide Ipi_blocks or Htt/Htp/Hpp")
    }
    blocks <- geoess_Ipi_hier(Htt, Htp, Hpp, return_components = TRUE)
    Ipi_eff <- blocks$Ipi_eff
    Ipi_between <- blocks$Ipi_between
    Ipi_borrowed <- blocks$Ipi_borrowed
  }

  if (psd != "none") {
    Ipi_eff <- .geoess_psd_adjust(Ipi_eff, psd)
  }

  I1_mat <- NULL
  if (is.function(I1)) {
    I1_mat <- I1(theta_ref, ...)
  } else if (is.list(I1) && !is.null(I1$type) && I1$type == "hessian") {
    if (is.null(I1$loglik_fn) || is.null(I1$n_obs)) {
      stop("I1 list must include loglik_fn and n_obs")
    }
    I1_mat <- geoess_I1_hessian(theta_ref, I1$loglik_fn, I1$n_obs, ...)
  } else {
    I1_mat <- I1
  }

  I1_mat <- .geoess_as_matrix(I1_mat, "I1")
  .geoess_check_square(I1_mat, "I1")

  if (method == "dir" && is.null(v)) stop("v is required for method = 'dir'")

  if (method == "trace") {
    ess_total <- geoess_trace(I1_mat, Ipi_eff)
    ess_between <- geoess_trace(I1_mat, Ipi_between)
    ess_borrowed <- geoess_trace(I1_mat, Ipi_borrowed)
  } else if (method == "dir") {
    ess_total <- geoess_dir(I1_mat, Ipi_eff, v)
    ess_between <- geoess_dir(I1_mat, Ipi_between, v)
    ess_borrowed <- geoess_dir(I1_mat, Ipi_borrowed, v)
  } else {
    ess_total <- geoess_eig(I1_mat, Ipi_eff)
    ess_between <- geoess_eig(I1_mat, Ipi_between)
    ess_borrowed <- geoess_eig(I1_mat, Ipi_borrowed)
  }

  out <- list(
    ess_total = ess_total,
    ess_between = ess_between,
    ess_borrowed = ess_borrowed,
    method = method,
    details = list(
      Ipi_eff = Ipi_eff,
      Ipi_between = Ipi_between,
      Ipi_borrowed = Ipi_borrowed,
      trace_additivity_check = if (method == "trace") {
        ess_total - (ess_between + ess_borrowed)
      } else {
        NA_real_
      }
    )
  )
  out
}

#' Hierarchical ESS from Stan Fits
#'
#' @param fit Stan fit object.
#' @param pars_theta Character vector of theta parameters.
#' @param pars_phi Character vector of phi parameters.
#' @param n_obs Integer number of observations.
#' @param I1_method Character. "hessian" or "score".
#' @param loglik_fn Log-likelihood function for Hessian route.
#' @param score_fn Score function for score route.
#' @param prior_source Character. "hierarchical_prior_fn" or "joint_draws_mix".
#' @param logprior_joint Joint log-prior function if using hierarchical_prior_fn.
#' @param Htt_fn,Htp_fn,Hpp_fn Optional block curvature functions.
#' @param mix_joint Optional pre-fit joint mixture.
#' @param ref Reference point choice for (theta, phi).
#' @param method ESS method.
#' @param ... Passed to internal computations.
#'
#' @return A list with hierarchical ESS results and reference points.
#' @export
geoess_hier_from_stan <- function(fit,
                                  pars_theta,
                                  pars_phi,
                                  n_obs,
                                  I1_method = c("hessian", "score"),
                                  loglik_fn = NULL,
                                  score_fn = NULL,
                                  prior_source = c("hierarchical_prior_fn", "joint_draws_mix"),
                                  logprior_joint = NULL,
                                  Htt_fn = NULL,
                                  Htp_fn = NULL,
                                  Hpp_fn = NULL,
                                  mix_joint = NULL,
                                  ref = c("posterior_mean", "posterior_mode"),
                                  method = c("trace", "eig", "dir"),
                                  ...) {
  I1_method <- match.arg(I1_method)
  prior_source <- match.arg(prior_source)
  ref <- match.arg(ref)
  method <- match.arg(method)

  ref_vals <- geoess_ref_hier(fit, pars_theta, pars_phi, ref = ref, ...)
  theta_ref <- ref_vals$theta_ref
  phi_ref <- ref_vals$phi_ref

  I1_mat <- NULL
  if (I1_method == "hessian") {
    if (is.null(loglik_fn)) stop("loglik_fn is required for Hessian I1")
    I1_mat <- geoess_I1_hessian(theta_ref, loglik_fn, n_obs, ...)
  } else {
    if (is.null(score_fn)) stop("score_fn is required for score I1")
    I1_mat <- geoess_I1_score(theta_ref, score_fn, n_obs, ...)
  }

  blocks <- NULL
  if (prior_source == "hierarchical_prior_fn") {
    blocks <- geoess_Ipi_hier_from_fn(
      theta_ref, phi_ref,
      Htt_fn = Htt_fn,
      Htp_fn = Htp_fn,
      Hpp_fn = Hpp_fn,
      logprior_joint = logprior_joint,
      ...
    )
  } else {
    if (is.null(mix_joint)) {
      draws_joint <- geoess_draws(fit, pars = c(pars_theta, pars_phi),
                                  type = "posterior", as = "matrix")
      mix_joint <- geoess_automixfit(draws_joint, family = "mvnorm",
                                     K_grid = 1:3, criterion = "BIC")
    }
    blocks <- geoess_Ipi_hier_from_mix(mix_joint, theta_ref, phi_ref)
  }

  res <- geoess_ess_hier(
    I1 = I1_mat,
    Ipi_blocks = blocks,
    theta_ref = theta_ref,
    phi_ref = phi_ref,
    method = method,
    ...
  )
  res$theta_ref <- theta_ref
  res$phi_ref <- phi_ref
  res$ref <- ref
  res
}
