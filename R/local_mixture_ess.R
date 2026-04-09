#' Local ESS for Mixture Priors with Component Attribution
#'
#' Computes local ESS for each mixture component at component-specific reference
#' points and reports optional attribution-scaled ESS.
#'
#' For a mixture prior \eqn{\pi(\theta) = \sum_k w_k \pi_k(\theta)}, this
#' function distinguishes:
#' \itemize{
#'   \item \strong{Component-local ESS} using component curvature
#'   \eqn{I_{\pi_k}(\theta^*_k)}.
#'   \item \strong{Mixture-local ESS} using full-mixture curvature
#'   \eqn{I_{\pi}(\theta^*_k)} (existing definition).
#' }
#'
#' Weight-scaled component reporting uses
#' \eqn{w_k \times \mathrm{ESS}_k(\theta_k^*)}. Responsibility-scaled reporting
#' uses \eqn{r_k(\theta_k^*) \times \mathrm{ESS}_k(\theta_k^*)}, where
#' \eqn{r_k(\theta) = \frac{w_k \pi_k(\theta)}{\sum_j w_j \pi_j(\theta)}}.
#' Hard-responsibility reporting uses \eqn{w_{k^*}\times \mathrm{ESS}_{k^*}},
#' where \eqn{k^*=\arg\max_k r_k(\theta^*)}. If ties occur, the hard
#' contribution is averaged over tied components:
#' \eqn{\frac{1}{|T|}\sum_{k\in T} w_k \times \mathrm{ESS}_k(\theta^*)},
#' with \eqn{T=\{k: r_k(\theta^*)=\max_j r_j(\theta^*)\}}.
#'
#' The returned \code{pointwise} table provides a unified strategy at each
#' evaluated reference point \eqn{\theta^*}:
#' \itemize{
#'   \item \eqn{\sum_k w_k \mathrm{ESS}_k(\theta^*)}
#'   \item \eqn{\sum_k r_k(\theta^*) \mathrm{ESS}_k(\theta^*)}
#'   \item \eqn{\sum_k h_k(\theta^*) \mathrm{ESS}_k(\theta^*)}, where
#'   \eqn{h_k(\theta^*) = w_k I(k\in T)/|T|} and \eqn{T} is the set of
#'   max-responsibility components
#'   \item \eqn{\mathrm{ESS}_{mix}(\theta^*)} from full-mixture curvature
#' }
#' so users can report attribution and full-mixture ESS together for all
#' reference points.
#'
#' @param mix Mixture prior object. Supports \code{geoess_mix} and
#'   \code{betaMix}. Custom objects may be used if they provide
#'   \code{weights}, and component/mixture log-density callbacks.
#' @param theta_star Optional reference point(s). If \code{NULL}, component-wise
#'   defaults are used (MVN mean/mode; Beta mode when interior, else mean).
#' @param component Character. \code{"all"} (default), \code{"max_resp"}, or
#'   \code{"given"}.
#' @param component_id Integer component index (or indices) when
#'   \code{component = "given"}.
#' @param I1 Unit Fisher information specification: matrix/scalar or function
#'   \code{I1(theta)}.
#' @param scale_by_weight Logical; if TRUE, include \eqn{w_k \times ESS_k}.
#' @param scale_by_responsibility Logical; if TRUE, include
#'   \eqn{r_k(\theta_k^*) \times ESS_k}.
#' @param ridge Numeric diagonal ridge added to \code{I1} for numerical
#'   stability.
#' @param use_cache Logical; if TRUE, use cached precision/Cholesky factors when
#'   available in \code{geoess_mix}.
#' @param symmetrize Logical; if TRUE, symmetrize curvature and Fisher matrices.
#' @param ... Passed to user-supplied \code{I1} or custom mixture callbacks.
#'
#' @return An object of class \code{geoess_local_mixture_result} with:
#' \itemize{
#'   \item \code{results}: data frame with one row per selected component and
#'   columns \code{k}, \code{weight}, \code{theta_star}, \code{ess_unscaled},
#'   \code{ess_scaled_weight}, \code{responsibility}, \code{ess_scaled_resp},
#'   \code{ess_mixture}.
#'   \item \code{pointwise}: one row per evaluated reference point containing
#'   unified local summaries from all components at that point:
#'   \code{ess_attr_weight}, \code{ess_attr_resp}, \code{ess_attr_hard_resp},
#'   and \code{ess_mixture}. This supports consistent reporting without
#'   point-specific special handling.
#'   \item \code{contributions}: long-form table of per-component contributions
#'   at each evaluated reference point.
#'   \item \code{settings}: list of call settings.
#' }
#'
#' @examples
#' mix <- structure(
#'   list(
#'     family = "mvnorm",
#'     weights = c(0.7, 0.3),
#'     mu = list(c(1), c(2)),
#'     Sigma = list(matrix(1, 1, 1), matrix(1, 1, 1))
#'   ),
#'   class = "geoess_mix"
#' )
#'
#' res <- geoess_local_mixture(
#'   mix = mix,
#'   theta_star = c(1, 2),
#'   I1 = matrix(1, 1, 1),
#'   scale_by_weight = TRUE
#' )
#' res$results
#' res$pointwise
#' @export
geoess_local_mixture <- function(mix,
                                 theta_star = NULL,
                                 component = c("all", "max_resp", "given"),
                                 component_id = NULL,
                                 I1,
                                 scale_by_weight = TRUE,
                                 scale_by_responsibility = FALSE,
                                 ridge = 0,
                                 use_cache = TRUE,
                                 symmetrize = TRUE,
                                 ...) {
  component <- match.arg(component)
  if (missing(I1)) stop("I1 is required")
  if (!is.logical(scale_by_weight) || length(scale_by_weight) != 1L) {
    stop("scale_by_weight must be TRUE/FALSE")
  }
  if (!is.logical(scale_by_responsibility) || length(scale_by_responsibility) != 1L) {
    stop("scale_by_responsibility must be TRUE/FALSE")
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0) {
    stop("ridge must be a nonnegative scalar")
  }

  w <- .geoess_mix_weights(mix)
  K <- length(w)

  if (component == "given") {
    if (is.null(component_id)) stop("component_id is required when component = 'given'")
    idx <- as.integer(component_id)
    if (length(idx) == 0L || any(!is.finite(idx)) || any(idx < 1L) || any(idx > K)) {
      stop("component_id must contain valid component indices")
    }
    idx <- unique(idx)
    theta_list <- .geoess_prepare_theta_star(theta_star, mix, idx = idx)
  } else if (component == "max_resp") {
    theta_one <- .geoess_prepare_theta_star_max_resp(theta_star, mix)
    resp <- .geoess_mix_responsibility(mix, theta_one, use_cache = use_cache, ...)
    idx <- which.max(resp)
    theta_list <- list(theta_one)
  } else {
    idx <- seq_len(K)
    theta_list <- .geoess_prepare_theta_star(theta_star, mix, idx = idx)
  }

  n_out <- length(idx)
  selected_attribution <- .geoess_local_mixture_selected_attribution(
    scale_by_weight = scale_by_weight,
    scale_by_responsibility = scale_by_responsibility
  )
  out <- data.frame(
    k = idx,
    weight = w[idx],
    ess_unscaled = rep(NA_real_, n_out),
    ess_scaled_weight = rep(NA_real_, n_out),
    responsibility = rep(NA_real_, n_out),
    ess_scaled_resp = rep(NA_real_, n_out),
    ess_mixture = rep(NA_real_, n_out)
  )
  out$theta_star <- I(vector("list", n_out))

  for (i in seq_len(n_out)) {
    k <- idx[i]
    theta_k <- theta_list[[i]]
    out$theta_star[[i]] <- theta_k

    I1_k <- .geoess_eval_I1(I1, theta_k, ridge = ridge, symmetrize = symmetrize, ...)
    Ipi_k <- .geoess_mix_component_Ipi(
      mix = mix,
      k = k,
      theta = theta_k,
      use_cache = use_cache,
      symmetrize = symmetrize,
      ...
    )

    ev_k <- eigen(.geoess_symmetrize(Ipi_k), symmetric = TRUE, only.values = TRUE)$values
    if (any(ev_k < -1e-8)) {
      warning("Component ", k, " prior curvature is indefinite; ESS may be negative")
    }

    ess_unscaled <- geoess_trace(I1_k, Ipi_k, symmetrize = symmetrize)
    out$ess_unscaled[i] <- ess_unscaled

    if (isTRUE(scale_by_weight)) {
      out$ess_scaled_weight[i] <- out$weight[i] * ess_unscaled
    }

    resp_k <- .geoess_mix_responsibility(mix, theta_k, use_cache = use_cache, ...)[k]
    out$responsibility[i] <- resp_k
    if (isTRUE(scale_by_responsibility)) {
      out$ess_scaled_resp[i] <- resp_k * ess_unscaled
    }

    Ipi_mix <- .geoess_mix_Ipi_full(
      mix = mix,
      theta = theta_k,
      use_cache = use_cache,
      symmetrize = symmetrize,
      ...
    )
    out$ess_mixture[i] <- geoess_trace(I1_k, Ipi_mix, symmetrize = symmetrize)
  }

  pointwise_obj <- .geoess_mix_pointwise_summary(
    mix = mix,
    theta_list = theta_list,
    I1 = I1,
    ridge = ridge,
    use_cache = use_cache,
    symmetrize = symmetrize,
    ...
  )
  pointwise_obj$pointwise$selected_attribution <- selected_attribution
  pointwise_obj$pointwise$ess_attributed <- switch(
    selected_attribution,
    weight = pointwise_obj$pointwise$ess_attr_weight,
    responsibility = pointwise_obj$pointwise$ess_attr_resp,
    none = rep(NA_real_, nrow(pointwise_obj$pointwise)),
    rep(NA_real_, nrow(pointwise_obj$pointwise))
  )

  res <- list(
    results = out,
    pointwise = pointwise_obj$pointwise,
    contributions = pointwise_obj$contributions,
    settings = list(
      component = component,
      component_id = idx,
      scale_by_weight = scale_by_weight,
      scale_by_responsibility = scale_by_responsibility,
      selected_attribution = selected_attribution,
      ridge = ridge,
      use_cache = use_cache
    ),
    call = match.call()
  )
  class(res) <- "geoess_local_mixture_result"
  res
}

#' @export
print.geoess_local_mixture_result <- function(x, digits = 4, ...) {
  cat("geoess_local_mixture_result\n")
  tab <- x$results
  theta_txt <- vapply(
    tab$theta_star,
    function(v) paste(format(signif(as.numeric(v), digits), trim = TRUE), collapse = ", "),
    character(1)
  )
  tab_print <- tab
  tab_print$theta_star <- theta_txt
  print(tab_print, row.names = FALSE, ...)

  if (!is.null(x$pointwise)) {
    cat("\npointwise summary (all components combined at each reference point)\n")
    ptab <- x$pointwise
    ptab$theta_star <- vapply(
      ptab$theta_star,
      function(v) paste(format(signif(as.numeric(v), digits), trim = TRUE), collapse = ", "),
      character(1)
    )
    print(ptab, row.names = FALSE, ...)
  }
  invisible(x)
}

.geoess_mix_weights <- function(mix) {
  w <- NULL
  if (inherits(mix, "geoess_mix")) {
    w <- mix$weights
  } else if (inherits(mix, "betaMix")) {
    w <- mix["w", ]
  } else if (is.list(mix) && !is.null(mix$weights)) {
    w <- mix$weights
  }
  if (is.null(w)) stop("Unable to extract mixture weights from `mix`")
  w <- as.numeric(w)
  if (any(!is.finite(w)) || any(w < 0)) stop("Mixture weights must be finite and nonnegative")
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) stop("Mixture weights must sum to a positive value")
  w / sw
}

.geoess_mix_dim <- function(mix) {
  if (inherits(mix, "geoess_mix")) {
    if (is.null(mix$mu) || length(mix$mu) == 0L) stop("mix$mu is missing")
    return(length(as.numeric(mix$mu[[1]])))
  }
  if (inherits(mix, "betaMix")) return(1L)
  if (is.list(mix) && !is.null(mix$dim)) return(as.integer(mix$dim)[1])
  stop("Unable to determine parameter dimension for `mix`")
}

.geoess_mix_default_theta_star <- function(mix, k) {
  if (inherits(mix, "geoess_mix")) {
    return(as.numeric(mix$mu[[k]]))
  }
  if (inherits(mix, "betaMix")) {
    a <- as.numeric(mix["a", k])
    b <- as.numeric(mix["b", k])
    if (is.finite(a) && is.finite(b) && a > 1 && b > 1) {
      return((a - 1) / (a + b - 2))
    }
    return(a / (a + b))
  }
  stop("theta_star must be supplied for this mixture type")
}

.geoess_mix_overall_mean <- function(mix) {
  w <- .geoess_mix_weights(mix)
  if (inherits(mix, "geoess_mix")) {
    mu <- do.call(rbind, lapply(mix$mu, as.numeric))
    return(as.numeric(crossprod(w, mu)))
  }
  if (inherits(mix, "betaMix")) {
    a <- as.numeric(mix["a", ])
    b <- as.numeric(mix["b", ])
    return(sum(w * (a / (a + b))))
  }
  if (is.list(mix) && !is.null(mix$mean)) return(as.numeric(mix$mean))
  stop("Unable to determine default overall mean for this mixture type")
}

.geoess_check_theta_dim <- function(theta, d) {
  theta <- as.numeric(theta)
  if (any(!is.finite(theta))) stop("theta_star contains non-finite values")
  if (length(theta) != d) stop("theta_star has incompatible dimension")
  theta
}

.geoess_prepare_theta_star <- function(theta_star, mix, idx = NULL) {
  w <- .geoess_mix_weights(mix)
  K <- length(w)
  d <- .geoess_mix_dim(mix)
  if (is.null(idx)) idx <- seq_len(K)

  if (is.null(theta_star)) {
    full <- lapply(seq_len(K), function(k) .geoess_mix_default_theta_star(mix, k))
    return(full[idx])
  }

  m <- length(idx)
  if (is.list(theta_star)) {
    if (length(theta_star) == m) {
      return(lapply(theta_star, .geoess_check_theta_dim, d = d))
    }
    if (length(theta_star) == K) {
      return(lapply(theta_star[idx], .geoess_check_theta_dim, d = d))
    }
    stop("theta_star list length must match selected components or all components")
  }

  if (is.matrix(theta_star)) {
    if (nrow(theta_star) == m) {
      return(lapply(seq_len(m), function(i) .geoess_check_theta_dim(theta_star[i, ], d = d)))
    }
    if (nrow(theta_star) == K) {
      return(lapply(idx, function(i) .geoess_check_theta_dim(theta_star[i, ], d = d)))
    }
    stop("theta_star matrix must have one row per selected/all component")
  }

  v <- as.numeric(theta_star)
  if (d == 1L && length(v) %in% c(m, K)) {
    if (length(v) == K) v <- v[idx]
    return(lapply(v, .geoess_check_theta_dim, d = d))
  }
  if (length(v) == d) {
    one <- .geoess_check_theta_dim(v, d = d)
    return(rep(list(one), m))
  }
  stop("Unable to parse theta_star; use vector, list, or matrix")
}

.geoess_prepare_theta_star_max_resp <- function(theta_star, mix) {
  d <- .geoess_mix_dim(mix)
  if (is.null(theta_star)) {
    return(.geoess_mix_overall_mean(mix))
  }
  if (is.list(theta_star)) {
    if (length(theta_star) < 1L) stop("theta_star list is empty")
    return(.geoess_check_theta_dim(theta_star[[1]], d = d))
  }
  if (is.matrix(theta_star)) {
    if (nrow(theta_star) < 1L) stop("theta_star matrix has zero rows")
    return(.geoess_check_theta_dim(theta_star[1, ], d = d))
  }
  .geoess_check_theta_dim(theta_star, d = d)
}

.geoess_mix_component_logdens <- function(mix, k, theta, use_cache = TRUE, ...) {
  theta <- as.numeric(theta)

  if (inherits(mix, "geoess_mix")) {
    mu_k <- as.numeric(mix$mu[[k]])
    Sig_k <- mix$Sigma[[k]]
    chol_k <- NULL
    if (isTRUE(use_cache) && !is.null(mix$chol) && length(mix$chol) >= k) {
      chol_k <- mix$chol[[k]]
    }
    return(.geoess_log_mvn(theta, mu_k, Sig_k, chol = chol_k))
  }

  if (inherits(mix, "betaMix")) {
    if (length(theta) != 1L) stop("betaMix component log density expects scalar theta")
    a <- as.numeric(mix["a", k])
    b <- as.numeric(mix["b", k])
    return(stats::dbeta(theta, a, b, log = TRUE))
  }

  if (is.list(mix) && is.function(mix$component_logdens)) {
    return(as.numeric(mix$component_logdens(k, theta, ...)))
  }

  stop("Component log-density is unavailable for this mixture object")
}

.geoess_mix_logdens_generic <- function(mix, theta, use_cache = TRUE, ...) {
  if (inherits(mix, "geoess_mix") && !is.null(mix$family) && mix$family == "mvnorm") {
    return(geoess_mix_logdens(mix, theta, use_cache = use_cache))
  }
  if (is.list(mix) && is.function(mix$logdens)) {
    return(as.numeric(mix$logdens(theta, ...)))
  }

  w <- .geoess_mix_weights(mix)
  log_comp <- vapply(
    seq_along(w),
    function(k) log(w[k]) + .geoess_mix_component_logdens(mix, k, theta, use_cache = use_cache, ...),
    numeric(1)
  )
  .geoess_logsumexp(log_comp)
}

.geoess_mix_responsibility <- function(mix, theta, use_cache = TRUE, ...) {
  w <- .geoess_mix_weights(mix)
  log_comp <- vapply(
    seq_along(w),
    function(k) log(w[k]) + .geoess_mix_component_logdens(mix, k, theta, use_cache = use_cache, ...),
    numeric(1)
  )
  .geoess_softmax_log(log_comp)
}

.geoess_mix_component_Ipi <- function(mix, k, theta, use_cache = TRUE, symmetrize = TRUE, ...) {
  if (inherits(mix, "geoess_mix")) {
    prec_k <- NULL
    if (isTRUE(use_cache) && !is.null(mix$prec) && length(mix$prec) >= k && !is.null(mix$prec[[k]])) {
      prec_k <- mix$prec[[k]]
    } else if (isTRUE(use_cache) && !is.null(mix$chol) && length(mix$chol) >= k && !is.null(mix$chol[[k]])) {
      prec_k <- .geoess_prec_from_chol(mix$chol[[k]])
    } else {
      prec_k <- .geoess_solve(mix$Sigma[[k]])
    }
    if (symmetrize) prec_k <- .geoess_symmetrize(prec_k)
    return(prec_k)
  }

  if (inherits(mix, "betaMix")) {
    a <- as.numeric(mix["a", k])
    b <- as.numeric(mix["b", k])
    out <- Ipi_beta(theta, a, b)
    if (symmetrize) out <- .geoess_symmetrize(out)
    return(out)
  }

  f <- function(z) .geoess_mix_component_logdens(mix, k, z, use_cache = use_cache, ...)
  H <- numDeriv::hessian(f, as.numeric(theta))
  Ipi <- -H
  if (symmetrize) Ipi <- .geoess_symmetrize(Ipi)
  Ipi
}

.geoess_mix_Ipi_full <- function(mix, theta, use_cache = TRUE, symmetrize = TRUE, ...) {
  if (inherits(mix, "geoess_mix") && !is.null(mix$family) && mix$family == "mvnorm") {
    return(geoess_Ipi_from_mix(
      mix = mix,
      theta = theta,
      use_cache = use_cache,
      symmetrize = symmetrize,
      psd = "none"
    ))
  }

  f <- function(z) .geoess_mix_logdens_generic(mix, z, use_cache = use_cache, ...)
  H <- numDeriv::hessian(f, as.numeric(theta))
  Ipi <- -H
  if (symmetrize) Ipi <- .geoess_symmetrize(Ipi)
  Ipi
}

.geoess_eval_I1 <- function(I1, theta, ridge = 0, symmetrize = TRUE, ...) {
  I1_mat <- if (is.function(I1)) I1(theta, ...) else I1
  I1_mat <- .geoess_as_matrix(I1_mat, "I1")
  .geoess_check_square(I1_mat, "I1")
  if (symmetrize) I1_mat <- .geoess_symmetrize(I1_mat)
  if (ridge > 0) {
    I1_mat <- I1_mat + diag(ridge, nrow(I1_mat))
  }
  I1_mat
}

.geoess_mix_attribution_scales <- function(mix,
                                           mixture_attribution = c("none", "weight", "responsibility", "hard_resp"),
                                           theta_star = NULL,
                                           component_idx = NULL,
                                           use_cache = TRUE,
                                           ...) {
  mixture_attribution <- match.arg(mixture_attribution)
  w <- .geoess_mix_weights(mix)
  K <- length(w)
  if (is.null(component_idx)) component_idx <- seq_len(K)
  component_idx <- as.integer(component_idx)

  if (mixture_attribution == "none") {
    return(rep(1, length(component_idx)))
  }
  if (mixture_attribution == "weight") {
    return(w[component_idx])
  }

  theta_list <- .geoess_prepare_theta_star(theta_star, mix, idx = component_idx)
  scales <- numeric(length(component_idx))
  for (i in seq_along(component_idx)) {
    k <- component_idx[i]
    r <- .geoess_mix_responsibility(mix, theta_list[[i]], use_cache = use_cache, ...)
    if (mixture_attribution == "responsibility") {
      scales[i] <- r[k]
    } else {
      hard <- .geoess_hard_resp_scales(r, w)
      scales[i] <- hard[k]
    }
  }
  scales
}

.geoess_hard_resp_scales <- function(resp, w, tol = 1e-10) {
  resp <- as.numeric(resp)
  w <- as.numeric(w)
  if (length(resp) != length(w)) stop("resp and w must have same length")
  max_resp <- max(resp)
  ties <- which(abs(resp - max_resp) <= tol)
  out <- numeric(length(resp))
  out[ties] <- w[ties] / length(ties)
  out
}

.geoess_local_mixture_selected_attribution <- function(scale_by_weight, scale_by_responsibility) {
  if (isTRUE(scale_by_responsibility)) return("responsibility")
  if (isTRUE(scale_by_weight)) return("weight")
  "none"
}

.geoess_mix_pointwise_summary <- function(mix,
                                          theta_list,
                                          I1,
                                          ridge = 0,
                                          use_cache = TRUE,
                                          symmetrize = TRUE,
                                          ...) {
  w <- .geoess_mix_weights(mix)
  K <- length(w)
  n_ref <- length(theta_list)

  pointwise <- data.frame(
    ref_id = seq_len(n_ref),
    k_max_resp = integer(n_ref),
    weight_max_resp = numeric(n_ref),
    responsibility_max_resp = numeric(n_ref),
    ess_attr_weight = numeric(n_ref),
    ess_attr_resp = numeric(n_ref),
    ess_attr_hard_resp = numeric(n_ref),
    ess_mixture = numeric(n_ref)
  )
  pointwise$theta_star <- I(vector("list", n_ref))

  contributions <- vector("list", n_ref)

  for (i in seq_len(n_ref)) {
    theta_i <- as.numeric(theta_list[[i]])
    pointwise$theta_star[[i]] <- theta_i

    I1_i <- .geoess_eval_I1(I1, theta_i, ridge = ridge, symmetrize = symmetrize, ...)
    ess_components <- numeric(K)

    for (k in seq_len(K)) {
      Ipi_k <- .geoess_mix_component_Ipi(
        mix = mix,
        k = k,
        theta = theta_i,
        use_cache = use_cache,
        symmetrize = symmetrize,
        ...
      )
      ess_components[k] <- geoess_trace(I1_i, Ipi_k, symmetrize = symmetrize)
    }

    resp <- .geoess_mix_responsibility(mix, theta_i, use_cache = use_cache, ...)
    hard <- .geoess_hard_resp_scales(resp, w)
    idx_hard <- which(hard > 0)
    if (!length(idx_hard)) idx_hard <- which.max(resp)
    k_hard <- idx_hard[1]

    contrib_weight <- w * ess_components
    contrib_resp <- resp * ess_components
    contrib_hard <- hard * ess_components

    Ipi_mix <- .geoess_mix_Ipi_full(
      mix = mix,
      theta = theta_i,
      use_cache = use_cache,
      symmetrize = symmetrize,
      ...
    )
    ess_mix <- geoess_trace(I1_i, Ipi_mix, symmetrize = symmetrize)

    pointwise$k_max_resp[i] <- k_hard
    pointwise$weight_max_resp[i] <- mean(w[idx_hard])
    pointwise$responsibility_max_resp[i] <- resp[k_hard]
    pointwise$ess_attr_weight[i] <- sum(contrib_weight)
    pointwise$ess_attr_resp[i] <- sum(contrib_resp)
    pointwise$ess_attr_hard_resp[i] <- sum(contrib_hard)
    pointwise$ess_mixture[i] <- ess_mix

    contributions[[i]] <- data.frame(
      ref_id = rep.int(i, K),
      k = seq_len(K),
      weight = w,
      responsibility = resp,
      ess_unscaled = ess_components,
      ess_contrib_weight = contrib_weight,
      ess_contrib_resp = contrib_resp,
      ess_contrib_hard_resp = contrib_hard
    )
  }

  contributions_df <- do.call(rbind, contributions)
  rownames(contributions_df) <- NULL

  list(pointwise = pointwise, contributions = contributions_df)
}
