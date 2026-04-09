#' Global Hierarchical ESS (Level-Specific)
#'
#' Computes global ESS for hyperparameters \eqn{\phi} and conditional group
#' parameters \eqn{\theta_i \mid \phi} by averaging local ESS over a weighting
#' distribution \eqn{Q(\Theta)}.
#'
#' @details
#' The hyperparameter ESS uses the local block ESS for \eqn{\phi} with a
#' per-group unit Fisher information. The conditional group ESS uses the local
#' conditional ESS for \eqn{\theta_i \mid \phi} obtained by Schur complements
#' on \eqn{I_1} and \eqn{I_\pi}. Global ESS is the expectation of these local
#' ESS values under the chosen weighting distribution.
#'
#' @param Q Character. "prior", "posterior", or "custom".
#' @param fit Optional Stan fit object for \code{Q = "prior"} or \code{Q = "posterior"}.
#' @param draws Optional draws matrix/data.frame of \eqn{\Theta} for \code{Q = "custom"}.
#' @param weights Optional nonnegative weights aligned with \code{draws}.
#' @param pars_theta Character vector of parameter names for \eqn{\theta}.
#' @param pars_phi Character vector of parameter names for \eqn{\phi}.
#' @param group_index Mapping of theta entries to groups. List of index vectors,
#'   or a vector of group labels of length \code{length(pars_theta)}.
#' @param Ipi_provider Character or function specifying how to obtain prior curvature.
#' @param Ipi_fn Function returning prior curvature blocks or full matrix.
#' @param Ipi_blocks Optional precomputed curvature blocks or arrays.
#' @param I1_provider Character or function specifying how to obtain unit Fisher info.
#' @param I1_theta_fn Function returning per-observation Fisher info for \eqn{\theta_i}.
#' @param I1_phi_fn Function returning per-group Fisher info for \eqn{\phi}.
#' @param I1_blocks Optional precomputed Fisher info blocks or arrays.
#' @param unit_phi Character. "group" only (default).
#' @param return_by_group Logical. If TRUE, return ESS per group.
#' @param psd Character. PSD handling for \eqn{I_\pi}.
#' @param ridge Numeric ridge added to diagonals for stability.
#' @param nsim Integer. Number of draws if sampling is needed.
#' @param seed Integer seed.
#' @param ... Passed to user-supplied functions.
#'
#' @return An object of class \code{geoess_hier_global}.
#' @examples
#' K <- 4
#' sigma <- 2
#' tau <- 0.5
#' s0 <- 1
#' draws <- cbind(matrix(rnorm(200 * K), ncol = K), rnorm(200))
#' pars_theta <- paste0("theta[", seq_len(K), "]")
#' pars_phi <- "mu"
#' colnames(draws) <- c(pars_theta, pars_phi)
#'
#' Ipi_fn <- function(Theta, ...) {
#'   Htt <- diag(1 / tau^2, K)
#'   Htp <- matrix(-1 / tau^2, nrow = K, ncol = 1)
#'   Hpp <- matrix(K / tau^2 + 1 / s0^2, nrow = 1, ncol = 1)
#'   list(Ipi_thetatheta = Htt, Ipi_thetaphi = Htp, Ipi_phiphi = Hpp)
#' }
#' I1_theta_fn <- function(theta_i, phi, ...) matrix(1 / sigma^2, 1, 1)
#' I1_phi_fn <- function(phi, theta, ...) matrix(1 / tau^2, 1, 1)
#'
#' res <- geoess_hier_global(
#'   Q = "custom",
#'   draws = draws,
#'   pars_theta = pars_theta,
#'   pars_phi = pars_phi,
#'   group_index = seq_len(K),
#'   Ipi_provider = "analytic",
#'   Ipi_fn = Ipi_fn,
#'   I1_provider = "analytic",
#'   I1_theta_fn = I1_theta_fn,
#'   I1_phi_fn = I1_phi_fn
#' )
#' res$ess_phi_global$estimate
#' res$ess_theta_cond_global
#' @export
geoess_hier_global <- function(Q = c("prior", "posterior", "custom"),
                               fit = NULL,
                               draws = NULL,
                               weights = NULL,
                               pars_theta,
                               pars_phi,
                               group_index = NULL,
                               Ipi_provider = c("analytic", "mixfit", "precomputed"),
                               Ipi_fn = NULL,
                               Ipi_blocks = NULL,
                               I1_provider = c("analytic", "precomputed"),
                               I1_theta_fn = NULL,
                               I1_phi_fn = NULL,
                               I1_blocks = NULL,
                               unit_phi = c("group"),
                               return_by_group = TRUE,
                               psd = c("none", "clip", "nearPD"),
                               ridge = 1e-8,
                               nsim = 2000,
                               seed = 1,
                               ...) {
  Q <- match.arg(Q)
  psd <- match.arg(psd)
  unit_phi <- match.arg(unit_phi)
  if (unit_phi != "group") stop("unit_phi must be 'group'")

  if (missing(pars_theta) || missing(pars_phi)) {
    stop("pars_theta and pars_phi are required")
  }
  pars_theta <- as.character(pars_theta)
  pars_phi <- as.character(pars_phi)

  if (is.function(Ipi_provider)) {
    Ipi_fn <- Ipi_provider
    Ipi_provider <- "analytic"
  }
  if (is.function(I1_provider)) {
    I1_fn <- I1_provider
    I1_provider <- "analytic"
  } else {
    I1_fn <- NULL
  }

  Ipi_provider <- match.arg(Ipi_provider)
  I1_provider <- match.arg(I1_provider)

  if (Q == "custom") {
    if (is.null(draws)) stop("draws are required for Q = 'custom'")
    draws_mat <- .geoess_hier_as_draws_matrix(draws, pars_theta, pars_phi)
  } else {
    if (!is.null(draws)) {
      draws_mat <- .geoess_hier_as_draws_matrix(draws, pars_theta, pars_phi)
    } else if (!is.null(fit)) {
      draws_mat <- geoess_draws(fit, pars = c(pars_theta, pars_phi),
                                type = Q, as = "matrix")
    } else {
      stop("fit or draws are required for Q = 'prior' or 'posterior'")
    }
  }

  n <- nrow(draws_mat)
  if (n < 1L) stop("No draws available")

  w_norm <- .geoess_normalize_weights(weights, n)

  group_list <- .geoess_group_index(group_index, pars_theta)
  K <- length(group_list)
  if (K < 1L) stop("No groups defined by group_index")

  idx_theta <- seq_len(length(pars_theta))
  idx_phi <- length(pars_theta) + seq_len(length(pars_phi))

  Ipi_get <- .geoess_Ipi_block_getter(Ipi_provider, Ipi_fn, Ipi_blocks,
                                      draws_mat, psd, idx_theta, idx_phi, ...)
  I1_get <- .geoess_I1_block_getter(I1_provider, I1_fn, I1_blocks,
                                    I1_theta_fn, I1_phi_fn, group_list,
                                    idx_theta, idx_phi, ...)

  ess_phi <- rep(NA_real_, n)
  ess_theta <- matrix(NA_real_, nrow = n, ncol = K)

  for (s in seq_len(n)) {
    Theta <- draws_mat[s, ]
    ipi_blocks <- Ipi_get(s, Theta)
    i1_blocks <- I1_get(s, Theta)

    ess_phi[s] <- .geoess_local_block_ess(i1_blocks$I1_phi, ipi_blocks$Ipi_phi, ridge = ridge)

    for (g in seq_len(K)) {
      idx_g <- group_list[[g]]
      idx_g_global <- idx_theta[idx_g]
      if (!is.null(ipi_blocks$Ipi_full)) {
        Ipi_cond <- .geoess_schur(ipi_blocks$Ipi_full, idx_g_global, idx_phi)
      } else {
        Ipi_cond <- .geoess_schur_blocks(
          ipi_blocks$Ipi_theta, ipi_blocks$Ipi_thetaphi, ipi_blocks$Ipi_phi,
          idx_g, ridge = ridge
        )
      }

      if (!is.null(i1_blocks$I1_full)) {
        I1_cond <- .geoess_schur(i1_blocks$I1_full, idx_g_global, idx_phi)
      } else if (!is.null(i1_blocks$I1_thetaphi)) {
        I1_cond <- .geoess_schur_blocks(
          i1_blocks$I1_theta, i1_blocks$I1_thetaphi, i1_blocks$I1_phi,
          idx_g, ridge = ridge
        )
      } else {
        I1_cond <- i1_blocks$I1_theta_list[[g]]
      }

      ess_theta[s, g] <- .geoess_local_block_ess(I1_cond, Ipi_cond, ridge = ridge)
    }
  }

  res_phi <- .geoess_global_summary(ess_phi, w_norm)
  res_theta <- apply(ess_theta, 2, .geoess_global_summary, w = w_norm)

  ess_theta_mean <- vapply(res_theta, `[[`, numeric(1), "estimate")
  ess_theta_mcse <- vapply(res_theta, `[[`, numeric(1), "mcse")

  out <- list(
    ess_phi_global = res_phi,
    ess_theta_cond_global = if (return_by_group) ess_theta_mean else mean(ess_theta_mean),
    ess_theta_mcse = if (return_by_group) ess_theta_mcse else mean(ess_theta_mcse),
    summaries = list(
      phi = res_phi,
      theta = if (return_by_group) res_theta else list(estimate = mean(ess_theta_mean))
    ),
    metadata = list(
      Q = Q,
      n_draws = n,
      unit_phi = unit_phi,
      pars_theta = pars_theta,
      pars_phi = pars_phi,
      group_index = group_list
    ),
    call = match.call()
  )
  class(out) <- "geoess_hier_global"
  out
}

.geoess_hier_as_draws_matrix <- function(draws, pars_theta, pars_phi) {
  if (is.data.frame(draws)) draws <- as.matrix(draws)
  if (!is.matrix(draws)) stop("draws must be a matrix or data.frame")

  if (!is.null(pars_theta) && !is.null(pars_phi)) {
    if (is.null(colnames(draws))) {
      if (ncol(draws) != length(pars_theta) + length(pars_phi)) {
        stop("draws has no column names and dimension mismatch with pars_theta/pars_phi")
      }
      return(draws)
    }
    idx <- match(c(pars_theta, pars_phi), colnames(draws))
    if (anyNA(idx)) stop("Missing parameters in draws: ", paste(c(pars_theta, pars_phi)[is.na(idx)], collapse = ", "))
    draws <- draws[, idx, drop = FALSE]
  }

  draws
}

.geoess_group_index <- function(group_index, pars_theta) {
  p <- length(pars_theta)
  if (p < 1L) stop("pars_theta must be non-empty")
  if (is.null(group_index)) {
    return(as.list(seq_len(p)))
  }
  if (is.list(group_index)) {
    return(lapply(group_index, function(idx) as.integer(idx)))
  }
  if (length(group_index) != p) stop("group_index length must match pars_theta")
  grp <- as.character(group_index)
  split(seq_len(p), grp)
}

.geoess_Ipi_block_getter <- function(Ipi_provider, Ipi_fn, Ipi_blocks,
                                     draws, psd, idx_theta, idx_phi, ...) {
  if (Ipi_provider == "mixfit") {
    mix <- geoess_automixfit(draws, family = "mvnorm", K_grid = 1:3, criterion = "BIC")
    return(function(i, Theta) {
      Ipi_full <- geoess_Ipi_from_mix(mix, Theta, psd = psd)
      list(
        Ipi_full = Ipi_full,
        Ipi_theta = Ipi_full[idx_theta, idx_theta, drop = FALSE],
        Ipi_thetaphi = Ipi_full[idx_theta, idx_phi, drop = FALSE],
        Ipi_phi = Ipi_full[idx_phi, idx_phi, drop = FALSE]
      )
    })
  }

  if (Ipi_provider == "precomputed") {
    if (is.null(Ipi_blocks)) stop("Ipi_blocks required for precomputed provider")
    Ipi_full_get <- .geoess_block_getter(Ipi_blocks$Ipi_full)
    Ipi_theta_get <- .geoess_block_getter(Ipi_blocks$Ipi_thetatheta)
    Ipi_thetaphi_get <- .geoess_block_getter(Ipi_blocks$Ipi_thetaphi)
    Ipi_phi_get <- .geoess_block_getter(Ipi_blocks$Ipi_phiphi)
    return(function(i, Theta) {
      Ipi_full <- Ipi_full_get(i)
      if (!is.null(Ipi_full)) {
        return(list(
          Ipi_full = Ipi_full,
          Ipi_theta = Ipi_full[idx_theta, idx_theta, drop = FALSE],
          Ipi_thetaphi = Ipi_full[idx_theta, idx_phi, drop = FALSE],
          Ipi_phi = Ipi_full[idx_phi, idx_phi, drop = FALSE]
        ))
      }
      list(
        Ipi_full = NULL,
        Ipi_theta = Ipi_theta_get(i),
        Ipi_thetaphi = Ipi_thetaphi_get(i),
        Ipi_phi = Ipi_phi_get(i)
      )
    })
  }

  if (is.function(Ipi_fn)) {
    return(function(i, Theta) {
      res <- Ipi_fn(Theta, ...)
      if (is.matrix(res)) {
        return(list(
          Ipi_full = res,
          Ipi_theta = res[idx_theta, idx_theta, drop = FALSE],
          Ipi_thetaphi = res[idx_theta, idx_phi, drop = FALSE],
          Ipi_phi = res[idx_phi, idx_phi, drop = FALSE]
        ))
      }
      if (is.list(res)) {
        if (!is.null(res$Ipi_full)) {
          full <- res$Ipi_full
          return(list(
            Ipi_full = full,
            Ipi_theta = full[idx_theta, idx_theta, drop = FALSE],
            Ipi_thetaphi = full[idx_theta, idx_phi, drop = FALSE],
            Ipi_phi = full[idx_phi, idx_phi, drop = FALSE]
          ))
        }
        return(list(
          Ipi_full = NULL,
          Ipi_theta = res$Ipi_thetatheta,
          Ipi_thetaphi = res$Ipi_thetaphi,
          Ipi_phi = res$Ipi_phiphi
        ))
      }
      stop("Ipi_fn must return a matrix or list of blocks")
    })
  }

  stop("Ipi_provider requires Ipi_fn, Ipi_blocks, or mixfit")
}

.geoess_I1_block_getter <- function(I1_provider, I1_fn, I1_blocks,
                                    I1_theta_fn, I1_phi_fn, group_list,
                                    idx_theta, idx_phi, ...) {
  if (I1_provider == "precomputed") {
    if (is.null(I1_blocks)) stop("I1_blocks required for precomputed provider")
    I1_full_get <- .geoess_block_getter(I1_blocks$I1_full)
    I1_theta_get <- .geoess_block_getter(I1_blocks$I1_thetatheta)
    I1_thetaphi_get <- .geoess_block_getter(I1_blocks$I1_thetaphi)
    I1_phi_get <- .geoess_block_getter(I1_blocks$I1_phiphi)
    return(function(i, Theta) {
      I1_full <- I1_full_get(i)
      if (!is.null(I1_full)) {
        return(list(
          I1_full = I1_full,
          I1_theta = I1_full[idx_theta, idx_theta, drop = FALSE],
          I1_thetaphi = I1_full[idx_theta, idx_phi, drop = FALSE],
          I1_phi = I1_full[idx_phi, idx_phi, drop = FALSE],
          I1_theta_list = NULL
        ))
      }
      list(
        I1_full = NULL,
        I1_theta = I1_theta_get(i),
        I1_thetaphi = I1_thetaphi_get(i),
        I1_phi = I1_phi_get(i),
        I1_theta_list = NULL
      )
    })
  }

  if (is.function(I1_fn)) {
    return(function(i, Theta) {
      res <- I1_fn(Theta, ...)
      if (is.matrix(res)) {
        return(list(
          I1_full = res,
          I1_theta = res[idx_theta, idx_theta, drop = FALSE],
          I1_thetaphi = res[idx_theta, idx_phi, drop = FALSE],
          I1_phi = res[idx_phi, idx_phi, drop = FALSE],
          I1_theta_list = NULL
        ))
      }
      if (is.list(res)) {
        if (!is.null(res$I1_full)) {
          full <- res$I1_full
          return(list(
            I1_full = full,
            I1_theta = full[idx_theta, idx_theta, drop = FALSE],
            I1_thetaphi = full[idx_theta, idx_phi, drop = FALSE],
            I1_phi = full[idx_phi, idx_phi, drop = FALSE],
            I1_theta_list = NULL
          ))
        }
        return(list(
          I1_full = NULL,
          I1_theta = res$I1_thetatheta,
          I1_thetaphi = res$I1_thetaphi,
          I1_phi = res$I1_phiphi,
          I1_theta_list = NULL
        ))
      }
      stop("I1_fn must return a matrix or list of blocks")
    })
  }

  if (is.null(I1_theta_fn) || is.null(I1_phi_fn)) {
    stop("I1_theta_fn and I1_phi_fn are required for analytic I1 provider")
  }

  return(function(i, Theta) {
    theta <- Theta[idx_theta]
    phi <- Theta[idx_phi]
    I1_phi <- .geoess_call_phi_theta(I1_phi_fn, phi, theta, ...)
    I1_phi <- .geoess_as_matrix(I1_phi, "I1_phi")

    I1_theta_list <- lapply(group_list, function(idx) {
      th_i <- theta[idx]
      val <- .geoess_call_theta_phi(I1_theta_fn, th_i, phi, ...)
      .geoess_as_matrix(val, "I1_theta")
    })

    list(
      I1_full = NULL,
      I1_theta = NULL,
      I1_thetaphi = NULL,
      I1_phi = I1_phi,
      I1_theta_list = I1_theta_list
    )
  })
}

.geoess_block_getter <- function(block) {
  if (is.null(block)) return(function(i) NULL)
  if (is.matrix(block) || inherits(block, "Matrix")) {
    block <- .geoess_as_matrix(block, "block")
    return(function(i) block)
  }
  if (is.list(block)) {
    return(function(i) .geoess_as_matrix(block[[i]], "block"))
  }
  if (is.array(block) && length(dim(block)) == 3L) {
    return(function(i) .geoess_as_matrix(block[, , i], "block"))
  }
  stop("Unsupported block type")
}

.geoess_schur_blocks <- function(M_theta, M_thetaphi, M_phi, idx_g, ridge = 1e-8) {
  if (is.null(M_theta) || is.null(M_thetaphi) || is.null(M_phi)) {
    stop("Block matrices are required for Schur complement")
  }
  M_aa <- M_theta[idx_g, idx_g, drop = FALSE]
  M_ab <- M_thetaphi[idx_g, , drop = FALSE]
  M_bb <- M_phi
  if (ridge > 0) {
    M_bb <- M_bb + diag(ridge, nrow(M_bb))
  }
  M_aa - M_ab %*% .geoess_solve(M_bb, t(M_ab))
}

.geoess_local_block_ess <- function(I1A, IpiA, ridge = 1e-8) {
  I1A <- .geoess_as_matrix(I1A, "I1A")
  IpiA <- .geoess_as_matrix(IpiA, "IpiA")
  I1A <- .geoess_symmetrize(I1A)
  IpiA <- .geoess_symmetrize(IpiA)
  if (ridge > 0) {
    I1A <- I1A + diag(ridge, nrow(I1A))
  }
  W <- .geoess_solve(I1A, IpiA)
  .geoess_trace(W) / nrow(I1A)
}

.geoess_global_summary <- function(x, w) {
  valid <- is.finite(x)
  if (!all(valid)) {
    x <- x[valid]
    w <- w[valid]
    w <- w / sum(w)
  }
  est <- .geoess_weighted_mean(x, w)
  mcse <- sqrt(.geoess_weighted_var(x, w)) / sqrt(.geoess_is_ess(w))
  quant_probs <- c(0.05, 0.5, 0.95)
  quant <- .geoess_weighted_quantile(x, w, probs = quant_probs)
  names(quant) <- paste0("q", quant_probs)
  list(estimate = est, mcse = mcse, quantiles = quant)
}
