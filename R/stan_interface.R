#' Extract Draws from Stan Fits
#'
#' @param fit A Stan fit object (rstan or cmdstanr), or a draws matrix/array.
#' @param pars Character vector of parameter names to extract.
#' @param type Character. "posterior" or "prior" (for metadata only).
#' @param thin Integer thinning interval.
#' @param as Character. "matrix" or "array".
#'
#' @return A matrix or array of draws.
#' @export
geoess_draws <- function(fit, pars = NULL, type = c("posterior", "prior"),
                         thin = 1, as = c("matrix", "array")) {
  type <- match.arg(type)
  as <- match.arg(as)
  if (!is.numeric(thin) || length(thin) != 1L || thin <= 0) {
    stop("thin must be a positive integer")
  }

  match_pars <- function(names_vec, p) {
    if (is.null(names_vec)) return(integer(0))
    exact <- which(names_vec == p)
    if (length(exact) > 0) return(exact)
    prefix <- which(startsWith(names_vec, paste0(p, "[")) |
                      startsWith(names_vec, paste0(p, ".")))
    if (length(prefix) > 0) return(prefix)
    integer(0)
  }

  if (inherits(fit, "stanfit")) {
    if (!requireNamespace("rstan", quietly = TRUE)) {
      stop("rstan is required for stanfit objects")
    }
    if (as == "matrix") {
      draws <- tryCatch(
        rstan::extract(fit, pars = pars, permuted = FALSE, inc_warmup = FALSE),
        error = function(e) NULL
      )
      if (is.null(draws) || isS4(draws) || !is.array(draws) || !is.numeric(draws)) {
        draws <- tryCatch(
          rstan::extract(fit, pars = pars, permuted = TRUE, inc_warmup = FALSE),
          error = function(e) NULL
        )
      }
      if (is.null(draws)) {
        draws <- base::as.matrix(fit)
      }
    } else {
      draws <- tryCatch(
        rstan::extract(fit, pars = pars, permuted = FALSE, inc_warmup = FALSE),
        error = function(e) NULL
      )
      if (is.null(draws) || isS4(draws) || !is.array(draws)) {
        draws <- tryCatch(
          rstan::extract(fit, pars = pars, permuted = TRUE, inc_warmup = FALSE),
          error = function(e) NULL
        )
      }
    }
  } else if (inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanVB")) {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      stop("cmdstanr is required for CmdStan* objects")
    }
    format <- if (as == "matrix") "draws_matrix" else "draws_array"
    draws <- fit$draws(variables = pars, format = format)
  } else if (is.matrix(fit) || is.array(fit)) {
    draws <- fit
  } else {
    stop("Unsupported fit type")
  }

  if (is.list(draws) && !is.array(draws) && !is.matrix(draws)) {
    draws <- .geoess_list_to_matrix(draws)
  }

  if (as == "matrix" && is.array(draws) && !is.matrix(draws)) {
    dims <- dim(draws)
    dimn <- dimnames(draws)
    if (!is.numeric(draws)) stop("Draws array must be numeric")
    draws <- matrix(draws, nrow = dims[1] * dims[2], ncol = dims[3])
    if (!is.null(dimn) && length(dimn) >= 3L) {
      colnames(draws) <- dimn[[3]]
    }
  }

  if (as == "array" && is.matrix(draws)) {
    draws <- array(draws, dim = c(nrow(draws), 1, ncol(draws)))
    if (!is.null(colnames(draws))) {
      dimnames(draws) <- list(NULL, NULL, colnames(draws))
    }
  }

  if (as == "matrix" && is.matrix(draws) && ncol(draws) == 0 && !is.null(pars)) {
    draws_all <- NULL
    if (inherits(fit, "stanfit")) {
      draws_all <- tryCatch(base::as.matrix(fit), error = function(e) NULL)
    } else if (inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanVB")) {
      draws_all <- tryCatch(fit$draws(variables = NULL, format = "draws_matrix"),
                            error = function(e) NULL)
    }
    if (!is.null(draws_all) && is.matrix(draws_all) && ncol(draws_all) > 0) {
      draws <- draws_all
    }
  }

  if (is.matrix(draws)) {
    if (!is.null(pars)) {
      if (is.character(pars)) {
        cn <- colnames(draws)
        if (is.null(cn)) stop("Draws matrix has no column names")
        idx <- integer(0)
        for (p in pars) {
          idx_p <- match_pars(cn, p)
          if (length(idx_p) > 0) {
            idx <- c(idx, idx_p)
          } else {
            any_match <- grep(p, cn, fixed = TRUE)
            if (length(any_match) > 0) {
              idx <- c(idx, any_match)
            } else {
              if (length(pars) == 1L) {
                warning("Parameter names not found in draws; returning all columns")
                idx <- seq_along(cn)
                break
              }
              stop("Missing parameters in draws: ", p)
            }
          }
        }
        draws <- draws[, idx, drop = FALSE]
      } else {
        draws <- draws[, pars, drop = FALSE]
      }
    }

    if (ncol(draws) == 0 && inherits(fit, "stanfit")) {
      draws_list <- tryCatch(
        rstan::extract(fit, pars = pars, permuted = TRUE, inc_warmup = FALSE),
        error = function(e) NULL
      )
      if (!is.null(draws_list)) {
        draws <- .geoess_list_to_matrix(draws_list)
      }
    }

    if (ncol(draws) == 0) {
      stop("No draws found for pars; check parameter names and Stan fit")
    }

    draws <- draws[seq(1, nrow(draws), by = thin), , drop = FALSE]
    return(draws)
  }

  if (is.array(draws)) {
    if (!is.null(pars)) {
      dimn <- dimnames(draws)
      if (!is.null(dimn) && length(dimn) >= 3L) {
        par_names <- dimn[[3]]
        if (is.character(pars)) {
          idx <- integer(0)
          for (p in pars) {
            idx_p <- match_pars(par_names, p)
            if (length(idx_p) > 0) {
              idx <- c(idx, idx_p)
            } else {
              any_match <- grep(p, par_names, fixed = TRUE)
              if (length(any_match) > 0) {
                idx <- c(idx, any_match)
              } else {
                if (length(pars) == 1L) {
                  warning("Parameter names not found in draws array; returning all elements")
                  idx <- seq_len(length(par_names))
                  break
                }
                stop("Missing parameters in draws array: ", p)
              }
            }
          }
        } else {
          idx <- pars
        }
        draws <- draws[, , idx, drop = FALSE]
      }
    }
    draws <- draws[seq(1, dim(draws)[1], by = thin), , , drop = FALSE]
    return(draws)
  }

  stop("Unexpected draws format")
}

.geoess_param_index_names <- function(base, dims) {
  if (length(dims) == 0L) return(base)
  idx <- expand.grid(lapply(dims, seq_len))
  idx_str <- apply(idx, 1L, paste, collapse = ",")
  paste0(base, "[", idx_str, "]")
}

.geoess_list_to_matrix <- function(draws) {
  if (is.numeric(draws) && !is.list(draws)) {
    return(matrix(draws, ncol = 1L))
  }
  if (!is.list(draws)) stop("draws must be a list or numeric vector")

  param_names <- names(draws)
  if (is.null(param_names) || any(param_names == "")) {
    param_names <- paste0("V", seq_along(draws))
  }

  mats <- vector("list", length(draws))
  col_names <- character(0)
  n_draws <- NULL

  for (i in seq_along(draws)) {
    x <- draws[[i]]
    if (!is.numeric(x)) stop("draws list elements must be numeric")

    d <- dim(x)
    if (is.null(d)) {
      mat <- matrix(x, ncol = 1L)
      cols <- param_names[i]
    } else {
      n_i <- d[1]
      n_cols <- prod(d[-1])
      mat <- matrix(x, nrow = n_i, ncol = n_cols)
      cols <- .geoess_param_index_names(param_names[i], d[-1])
    }

    if (is.null(n_draws)) {
      n_draws <- nrow(mat)
    } else if (nrow(mat) != n_draws) {
      stop("Inconsistent number of draws across parameters")
    }

    mats[[i]] <- mat
    col_names <- c(col_names, cols)
  }

  out <- do.call(cbind, mats)
  colnames(out) <- col_names
  out
}

#' Extract Log-Likelihood Draws
#'
#' Convenience wrapper for generated quantities named \code{log_lik}.
#'
#' @param fit A Stan fit object or draws matrix/array.
#' @param var Character name for log-likelihood variable (default "log_lik").
#' @param thin Integer thinning interval.
#' @param as Character. "matrix" or "array".
#'
#' @return A matrix or array of log-likelihood draws.
#' @export
geoess_loglik_draws <- function(fit, var = "log_lik", thin = 1,
                                as = c("matrix", "array")) {
  geoess_draws(fit, pars = var, type = "posterior", thin = thin, as = as)
}

.geoess_theta_ref_map <- function(fit, pars, ...) {
  if (inherits(fit, "stanfit")) {
    return(.geoess_theta_ref_map_rstan(fit, pars, ...))
  }
  if (inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanVB")) {
    return(.geoess_theta_ref_map_cmdstanr(fit, pars, ...))
  }
  NULL
}

.geoess_theta_ref_map_rstan <- function(fit, pars, ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) return(NULL)
  if (is.null(fit@stanmodel) || is.null(fit@data)) return(NULL)

  opt <- tryCatch(
    rstan::optimizing(fit@stanmodel, data = fit@data, ...),
    error = function(e) NULL
  )
  if (is.null(opt) || is.null(opt$par)) return(NULL)

  par_vec <- unlist(opt$par, use.names = TRUE)
  if (!is.null(pars)) {
    missing <- setdiff(pars, names(par_vec))
    if (length(missing) > 0) return(NULL)
    par_vec <- par_vec[pars]
  }
  as.numeric(par_vec)
}

.geoess_theta_ref_map_cmdstanr <- function(fit, pars, ...) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(NULL)
  model <- tryCatch(fit$cmdstan_model(), error = function(e) NULL)
  if (is.null(model)) return(NULL)

  data <- NULL
  meta <- tryCatch(fit$metadata(), error = function(e) NULL)
  if (!is.null(meta) && !is.null(meta$data)) data <- meta$data

  if (is.null(data)) return(NULL)

  opt <- tryCatch(model$optimize(data = data, ...), error = function(e) NULL)
  if (is.null(opt)) return(NULL)

  summ <- tryCatch(opt$summary(), error = function(e) NULL)
  if (is.null(summ)) return(NULL)

  if (!all(c("variable", "estimate") %in% names(summ))) {
    if (all(c("variable", "mean") %in% names(summ))) {
      summ$estimate <- summ$mean
    } else {
      return(NULL)
    }
  }

  par_vec <- summ$estimate
  names(par_vec) <- summ$variable
  if (!is.null(pars)) {
    missing <- setdiff(pars, names(par_vec))
    if (length(missing) > 0) return(NULL)
    par_vec <- par_vec[pars]
  }
  as.numeric(par_vec)
}
