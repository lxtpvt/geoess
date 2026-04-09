#' Geometric Effective Sample Size from Curvature
#'
#' Computes local, coordinate-free ESS from Fisher information and prior curvature.
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature (negative Hessian of log prior).
#' @param symmetrize Logical. If TRUE, symmetrize inputs before eigen calculations.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{ess}: scalar ESS, \eqn{\frac{1}{d}\,\mathrm{tr}(I_1^{-1} I_\pi)}.
#'   \item \code{eigen_ess}: eigenvalues of the whitened prior curvature.
#'   \item \code{dir_ess}: function(v) computing directional ESS for vector v.
#' }
#'
#' @examples
#' I1 <- diag(2)
#' Ipi <- diag(c(0.5, 0.25))
#' geoess_curvature(I1, Ipi)$ess
#' @export
geoess_curvature <- function(I1, Ipi, symmetrize = TRUE) {
  I1 <- .geoess_as_matrix(I1, "I1")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I1, "I1")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I1, Ipi, "I1", "Ipi")

  if (symmetrize) {
    I1 <- .geoess_symmetrize(I1)
    Ipi <- .geoess_symmetrize(Ipi)
  }

  ess_scalar <- geoess_trace(I1, Ipi, symmetrize = symmetrize)
  eigen_ess <- geoess_eigen(I1, Ipi, symmetrize = symmetrize)
  dir_ess <- function(v) geoess_directional(I1, Ipi, v)

  list(ess = ess_scalar, eigen_ess = eigen_ess, dir_ess = dir_ess)
}

#' Directional Effective Sample Size
#'
#' Computes ESS for a specific direction in parameter space.
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature, or a \code{geoess_mix} object for mixtures.
#' @param v Numeric vector. Direction in parameter space.
#'
#' @return Numeric scalar ESS in direction v.
#'
#' @examples
#' I1 <- diag(2)
#' Ipi <- diag(c(1, 3))
#' geoess_directional(I1, Ipi, c(1, 0))
#' @export
geoess_directional <- function(I1, Ipi, v) {
  I1 <- .geoess_as_matrix(I1, "I1")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I1, "I1")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I1, Ipi, "I1", "Ipi")

  v <- as.numeric(v)
  if (length(v) != nrow(I1)) stop("Length of v must match dimension of I1")
  if (any(!is.finite(v))) stop("v contains non-finite entries")

  numerator <- as.numeric(t(v) %*% Ipi %*% v)
  denominator <- as.numeric(t(v) %*% I1 %*% v)
  if (denominator <= 0) stop("Directional denominator is nonpositive")

  numerator / denominator
}

#' Eigen ESS Spectrum
#'
#' Computes the eigen-ESS spectrum using a whitened prior curvature.
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param symmetrize Logical. If TRUE, symmetrize inputs before eigen calculations.
#'
#' @return Numeric vector of eigenvalues.
#'
#' @examples
#' I1 <- diag(c(2, 3))
#' Ipi <- diag(c(1, 6))
#' geoess_eigen(I1, Ipi)
#' @export
geoess_eigen <- function(I1, Ipi, symmetrize = TRUE) {
  I1 <- .geoess_as_matrix(I1, "I1")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I1, "I1")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I1, Ipi, "I1", "Ipi")

  if (symmetrize) {
    I1 <- .geoess_symmetrize(I1)
    Ipi <- .geoess_symmetrize(Ipi)
  }

  I1_dense <- as.matrix(I1)
  Ipi_dense <- as.matrix(Ipi)
  chol_fit <- .geoess_try_chol(I1_dense)

  if (is.null(chol_fit$R)) {
    warning("I1 is not positive definite; using symmetrized I1^{-1} Ipi")
    W <- .geoess_solve(I1_dense, Ipi_dense)
    W <- .geoess_symmetrize(W)
  } else {
    if (chol_fit$jittered) warning("I1 was jittered for Cholesky stability")
    R <- chol_fit$R
    A <- backsolve(R, diag(nrow(I1_dense)))
    W <- A %*% Ipi_dense %*% t(A)
    W <- .geoess_symmetrize(W)
  }

  eigen(W, symmetric = TRUE, only.values = TRUE)$values
}

#' ESS Field over Contexts
#'
#' Applies ESS computations to a list of Fisher information matrices.
#'
#' @param I1_list List of Fisher information matrices.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param method Character. One of "trace", "eig", or "dir".
#' @param v Optional numeric vector for directional ESS when method = "dir".
#' @param symmetrize Logical. If TRUE, symmetrize inputs before eigen calculations.
#'
#' @return A list or numeric vector of ESS values by context.
#'
#' @examples
#' I1_list <- list(diag(2), diag(c(2, 3)))
#' Ipi <- diag(2)
#' geoess_field(I1_list, Ipi, method = "trace")
#' @export
geoess_field <- function(I1_list, Ipi, method = c("trace", "eig", "dir"), v = NULL,
                         symmetrize = TRUE) {
  method <- match.arg(method)
  if (!is.list(I1_list) || length(I1_list) == 0) {
    stop("I1_list must be a non-empty list of matrices")
  }

  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(Ipi, "Ipi")

  if (method == "dir" && is.null(v)) stop("v is required when method = 'dir'")

  res <- lapply(I1_list, function(I1) {
    if (method == "trace") {
      geoess_trace(I1, Ipi, symmetrize = symmetrize)
    } else if (method == "eig") {
      geoess_eigen(I1, Ipi, symmetrize = symmetrize)
    } else {
      geoess_directional(I1, Ipi, v)
    }
  })

  if (method %in% c("trace", "dir")) unlist(res) else res
}

#' Trace ESS (Scalar)
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param symmetrize Logical. If TRUE, symmetrize inputs before calculations.
#'
#' @return Numeric scalar ESS.
#' @export
geoess_trace <- function(I1, Ipi, symmetrize = TRUE) {
  I1 <- .geoess_as_matrix(I1, "I1")
  Ipi <- .geoess_as_matrix(Ipi, "Ipi")
  .geoess_check_square(I1, "I1")
  .geoess_check_square(Ipi, "Ipi")
  .geoess_check_same_dim(I1, Ipi, "I1", "Ipi")

  if (symmetrize) {
    I1 <- .geoess_symmetrize(I1)
    Ipi <- .geoess_symmetrize(Ipi)
  }

  d <- nrow(I1)
  W <- .geoess_solve(I1, Ipi)
  .geoess_trace(W) / d
}

geoess_ess_trace <- function(I1, Ipi, symmetrize = TRUE) {
  geoess_trace(I1, Ipi, symmetrize = symmetrize)
}

geoess_ess_dir <- function(I1, Ipi, v) {
  geoess_directional(I1, Ipi, v)
}

geoess_ess_eigs <- function(I1, Ipi, symmetrize = TRUE) {
  geoess_eigen(I1, Ipi, symmetrize = symmetrize)
}

#' Directional ESS (Alias)
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param v Numeric vector direction.
#'
#' @return Numeric scalar ESS for direction v.
#' @export
geoess_dir <- function(I1, Ipi, v) {
  geoess_directional(I1, Ipi, v)
}

#' Eigen ESS Spectrum (Alias)
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param symmetrize Logical. If TRUE, symmetrize inputs before eigen calculations.
#'
#' @return Numeric vector of eigenvalues.
#' @export
geoess_eig <- function(I1, Ipi, symmetrize = TRUE) {
  geoess_eigen(I1, Ipi, symmetrize = symmetrize)
}

#' ESS Dispatcher
#'
#' @param I1 Numeric matrix. Fisher information per observation.
#' @param Ipi Numeric matrix. Prior curvature.
#' @param method Character. One of "trace", "eig", "dir", or "hier".
#' @param v Numeric vector for directional ESS (method = "dir").
#' @param Htt Numeric matrix for target block (method = "hier").
#' @param Htp Numeric matrix for cross block (method = "hier").
#' @param Hpp Numeric matrix for nuisance block (method = "hier").
#' @param symmetrize Logical. If TRUE, symmetrize inputs before eigen calculations.
#' @param mixture_attribution Character. For mixture priors in local ESS:
#'   "none" (default), "weight", "responsibility", or "hard_resp". Only
#'   applied when
#'   \code{Ipi} is a \code{geoess_mix} object. For \code{"hard_resp"}, local
#'   ESS uses weighted hard assignment: \eqn{w_{k^*}\times ESS_{k^*}} for
#'   unique \eqn{k^*=\arg\max_k r_k(\theta^*)}; ties are averaged across
#'   max-responsibility components. Ignored when
#'   \code{mixture_curvature = "mix"}.
#' @param mixture_curvature Character. For mixture priors in local ESS:
#'   \code{"component"} (default) uses component curvature and preserves the
#'   existing weighted component ESS behavior; \code{"mix"} uses full-mixture
#'   curvature \eqn{-\nabla^2 \log \sum_k w_k \pi_k(\theta^*)} at
#'   \code{theta_star}.
#' @param theta_star Optional reference point(s) used when
#'   \code{mixture_attribution = "responsibility"},
#'   \code{mixture_attribution = "hard_resp"}, or
#'   \code{mixture_curvature = "mix"}. For \code{mixture_curvature = "mix"},
#'   \code{theta_star} must be a single numeric vector.
#'
#' @return ESS result, mixture ESS list when \code{Ipi} is a mixture, or
#'   Schur complement for method = "hier".
#' @export
geoess_ess <- function(I1 = NULL, Ipi = NULL,
                       method = c("trace", "eig", "dir", "hier"),
                       v = NULL, Htt = NULL, Htp = NULL, Hpp = NULL,
                       symmetrize = TRUE,
                       mixture_attribution = c("none", "weight", "responsibility", "hard_resp"),
                       mixture_curvature = c("component", "mix"),
                       theta_star = NULL) {
  method <- match.arg(method)
  mixture_attribution <- match.arg(mixture_attribution)
  mixture_curvature <- match.arg(mixture_curvature)

  if (!is.null(Ipi) && inherits(Ipi, "geoess_mix") && method %in% c("trace", "dir")) {
    mix <- Ipi
    if (is.null(mix$weights) || is.null(mix$Sigma)) {
      stop("mix must contain weights and Sigma")
    }

    I1 <- .geoess_as_matrix(I1, "I1")
    .geoess_check_square(I1, "I1")
    if (symmetrize) I1 <- .geoess_symmetrize(I1)

    w <- as.numeric(mix$weights)
    if (any(!is.finite(w)) || any(w < 0)) stop("mix$weights must be finite and nonnegative")
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) stop("mix$weights must sum to a positive value")
    w <- w / sw

    K <- length(w)
    d <- nrow(I1)
    if (length(mix$Sigma) != K) stop("mix$Sigma length must match mix$weights")

    if (method == "dir") {
      if (is.null(v)) stop("v is required for directional ESS")
      v <- as.numeric(v)
      if (length(v) != d) stop("Length of v must match dimension of I1")
      denom <- as.numeric(t(v) %*% I1 %*% v)
      if (denom <= 0) stop("Directional denominator is nonpositive")
    }

    if (mixture_curvature == "mix") {
      if (mixture_attribution != "none") {
        stop("mixture_attribution must be 'none' when mixture_curvature = 'mix'")
      }
      if (is.null(theta_star)) {
        stop("theta_star is required when mixture_curvature = 'mix'")
      }
      if (is.list(theta_star) || is.matrix(theta_star) || is.data.frame(theta_star)) {
        stop("theta_star must be a single numeric vector when mixture_curvature = 'mix'")
      }

      theta_one <- as.numeric(theta_star)
      if (any(!is.finite(theta_one))) {
        stop("theta_star contains non-finite values")
      }
      if (length(theta_one) != d) {
        stop("Length of theta_star must match dimension of I1")
      }

      Ipi_mix <- geoess_Ipi_from_mix(
        mix = mix,
        theta = theta_one,
        symmetrize = symmetrize,
        psd = "none"
      )

      ess_mix <- if (method == "trace") {
        geoess_trace(I1, Ipi_mix, symmetrize = symmetrize)
      } else {
        geoess_directional(I1, Ipi_mix, v)
      }

      return(list(
        ess = ess_mix,
        ess_mixture = ess_mix,
        Ipi_mixture = Ipi_mix,
        theta_star = theta_one,
        mixture_curvature = mixture_curvature,
        method = method
      ))
    }

    ess_components <- numeric(K)
    for (k in seq_len(K)) {
      Sig_k <- mix$Sigma[[k]]
      if (!is.matrix(Sig_k) || any(dim(Sig_k) != c(d, d))) {
        stop("mix$Sigma[[k]] has incompatible dimensions")
      }
      prec_k <- NULL
      if (!is.null(mix$prec) && length(mix$prec) >= k && !is.null(mix$prec[[k]])) {
        prec_k <- mix$prec[[k]]
      } else if (!is.null(mix$chol) && length(mix$chol) >= k && !is.null(mix$chol[[k]])) {
        prec_k <- .geoess_prec_from_chol(mix$chol[[k]])
      } else {
        prec_k <- .geoess_solve(Sig_k)
      }

      if (symmetrize) prec_k <- .geoess_symmetrize(prec_k)

      if (method == "trace") {
        W <- .geoess_solve(I1, prec_k)
        ess_components[k] <- .geoess_trace(W) / d
      } else {
        num <- as.numeric(t(v) %*% prec_k %*% v)
        ess_components[k] <- num / denom
      }
    }

    ess_weighted <- sum(w * ess_components)
    attr_scales <- rep(1, K)
    ess_attributed <- ess_weighted
    ess_components_attributed <- ess_components

    if (mixture_attribution != "none") {
      attr_scales <- .geoess_mix_attribution_scales(
        mix = mix,
        mixture_attribution = mixture_attribution,
        theta_star = theta_star,
        component_idx = seq_len(K)
      )
      ess_components_attributed <- ess_components * attr_scales
      ess_attributed <- sum(ess_components_attributed)
    }

    return(list(
      ess = if (mixture_attribution == "none") ess_weighted else ess_attributed,
      ess_weighted = ess_weighted,
      ess_components = ess_components,
      ess_components_attributed = ess_components_attributed,
      attribution_scales = attr_scales,
      mixture_attribution = mixture_attribution,
      mixture_curvature = mixture_curvature,
      ess_attributed = ess_attributed,
      weights = w,
      method = method
    ))
  } else if (!is.null(Ipi) && inherits(Ipi, "geoess_mix")) {
    stop("method '", method, "' not supported for mixture prior")
  }

  if (method == "trace") {
    return(geoess_trace(I1, Ipi, symmetrize = symmetrize))
  }
  if (method == "eig") {
    return(geoess_eig(I1, Ipi, symmetrize = symmetrize))
  }
  if (method == "dir") {
    if (is.null(v)) stop("v is required for directional ESS")
    return(geoess_dir(I1, Ipi, v))
  }

  if (is.null(Htt) || is.null(Htp) || is.null(Hpp)) {
    stop("Htt, Htp, and Hpp are required for hierarchical ESS")
  }
  schur_theta_given_phi(Htt, Htp, Hpp, return_components = FALSE)
}
