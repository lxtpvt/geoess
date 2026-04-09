#' Plot 1D Parameter Draws
#'
#' Visualizes 1D parameter draws from a Stan fit or a numeric matrix.
#'
#' @param x Stan fit, numeric matrix/array, or numeric vector of draws.
#' @param pars Optional character vector of parameter names to extract.
#' @param type Character. "posterior" or "prior" (metadata only).
#' @param thin Integer thinning interval.
#' @param method Character. "hist", "density", or "both".
#' @param bins Integer number of histogram bins.
#' @param main Optional plot title(s).
#' @param col Histogram fill color.
#' @param border Histogram border color.
#' @param lwd Line width for density curves.
#' @param ... Passed to plotting functions.
#'
#' @return Invisible list of histogram/density objects by parameter.
#'
#' @examples
#' set.seed(123)
#' draws <- cbind(rnorm(200))
#' geoess_plot_draws_1d(draws, method = "both")
#' @export
geoess_plot_draws_1d <- function(x, pars = NULL, type = c("posterior", "prior"),
                                 thin = 1, method = c("both", "hist", "density"),
                                 bins = 30, main = NULL, col = "gray80",
                                 border = "white", lwd = 2, ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  if (!is.numeric(thin) || length(thin) != 1L || thin <= 0) {
    stop("thin must be a positive integer")
  }
  if (!is.numeric(bins) || length(bins) != 1L || bins <= 0) {
    stop("bins must be a positive integer")
  }

  if (is.numeric(x) && is.null(dim(x))) {
    if (!is.null(pars) && length(pars) != 1L) {
      stop("pars must be length 1 for vector draws")
    }
    draws <- matrix(x, ncol = 1L)
    colnames(draws) <- if (!is.null(pars)) as.character(pars) else "V1"
  } else {
    draws <- geoess_draws(x, pars = pars, type = type, thin = thin, as = "matrix")
  }

  if (!is.matrix(draws) || ncol(draws) == 0L) {
    stop("No draws available for plotting")
  }

  par_names <- colnames(draws)
  if (is.null(par_names)) par_names <- paste0("V", seq_len(ncol(draws)))

  n_par <- ncol(draws)
  if (n_par > 1L) {
    ncol_plot <- if (n_par <= 3L) n_par else 2L
    nrow_plot <- ceiling(n_par / ncol_plot)
    op <- par(mfrow = c(nrow_plot, ncol_plot), mar = c(4, 4, 2, 1))
    on.exit(par(op), add = TRUE)
  }

  res <- vector("list", n_par)
  for (i in seq_len(n_par)) {
    xi <- draws[, i]
    title_i <- if (is.null(main)) {
      par_names[i]
    } else {
      main[(i - 1L) %% length(main) + 1L]
    }
    if (method == "hist") {
      h <- hist(xi, breaks = bins, col = col, border = border,
                main = title_i, xlab = par_names[i], ...)
      res[[i]] <- list(hist = h)
    } else if (method == "density") {
      d <- stats::density(xi)
      plot(d, main = title_i, xlab = par_names[i], lwd = lwd, ...)
      res[[i]] <- list(density = d)
    } else {
      h <- hist(xi, breaks = bins, col = col, border = border, freq = FALSE,
                main = title_i, xlab = par_names[i], ...)
      d <- stats::density(xi)
      lines(d, lwd = lwd)
      res[[i]] <- list(hist = h, density = d)
    }
  }

  names(res) <- par_names
  invisible(res)
}

#' Plot 2D Parameter Draws
#'
#' Visualizes 2D parameter draws using scatter or density contours.
#'
#' @param x Stan fit or numeric matrix/array of draws.
#' @param pars Character vector of length 2 naming parameters to plot.
#' @param type Character. "posterior" or "prior" (metadata only).
#' @param thin Integer thinning interval.
#' @param method Character. "scatter", "smooth", or "contour".
#' @param n Integer grid size for contour density.
#' @param col Point color for scatter.
#' @param pch Point character for scatter.
#' @param cex Point size for scatter.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param ... Passed to plotting functions.
#'
#' @return Invisible list with x, y, and density (if computed).
#'
#' @examples
#' set.seed(123)
#' draws <- cbind(mu = rnorm(300), tau = rnorm(300, 1, 0.2))
#' geoess_plot_draws_2d(draws, pars = c("mu", "tau"), method = "scatter")
#' @export
geoess_plot_draws_2d <- function(x, pars = NULL, type = c("posterior", "prior"),
                                 thin = 1, method = c("scatter", "smooth", "contour"),
                                 n = 60, col = "gray30", pch = 16, cex = 0.6,
                                 xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  if (!is.numeric(thin) || length(thin) != 1L || thin <= 0) {
    stop("thin must be a positive integer")
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 10) {
    stop("n must be an integer > 10")
  }

  draws <- geoess_draws(x, pars = pars, type = type, thin = thin, as = "matrix")
  if (!is.matrix(draws) || ncol(draws) < 2L) {
    stop("draws must have at least two parameters for 2D plotting")
  }

  if (is.null(pars)) {
    idx <- 1:2
  } else if (is.character(pars)) {
    if (is.null(colnames(draws))) {
      stop("draws matrix has no column names; use numeric indices for pars")
    }
    idx <- match(pars, colnames(draws))
    if (any(is.na(idx))) stop("Could not find pars in draws")
  } else {
    idx <- as.integer(pars)
  }

  if (length(idx) != 2L) stop("pars must be length 2")

  xdat <- draws[, idx[1]]
  ydat <- draws[, idx[2]]

  if (is.null(xlab)) {
    xlab <- if (!is.null(colnames(draws))) colnames(draws)[idx[1]] else paste0("V", idx[1])
  }
  if (is.null(ylab)) {
    ylab <- if (!is.null(colnames(draws))) colnames(draws)[idx[2]] else paste0("V", idx[2])
  }

  out <- list(x = xdat, y = ydat, density = NULL)

  if (method == "scatter") {
    plot(xdat, ydat, xlab = xlab, ylab = ylab, col = col, pch = pch,
         cex = cex, ...)
  } else if (method == "smooth") {
    smoothScatter(xdat, ydat, xlab = xlab, ylab = ylab, ...)
  } else {
    kde <- .geoess_kde2d(xdat, ydat, n = n)
    contour(kde$x, kde$y, kde$z, xlab = xlab, ylab = ylab, ...)
    out$density <- kde
  }

  invisible(out)
}

.geoess_kde2d <- function(x, y, n = 60, lims = NULL, h = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) stop("x and y must be the same length")
  if (length(x) < 2L) stop("x and y must have at least two points")

  if (is.null(lims)) {
    lims <- c(range(x), range(y))
  }
  if (!is.numeric(lims) || length(lims) != 4L) {
    stop("lims must be a numeric vector of length 4")
  }

  if (is.null(h)) {
    h <- c(stats::bw.nrd(x), stats::bw.nrd(y))
  }
  if (any(!is.finite(h)) || any(h <= 0)) stop("invalid bandwidths")

  gx <- seq(lims[1], lims[2], length.out = n)
  gy <- seq(lims[3], lims[4], length.out = n)

  dx <- outer(gx, x, "-") / h[1]
  dy <- outer(gy, y, "-") / h[2]
  kx <- exp(-0.5 * dx^2)
  ky <- exp(-0.5 * dy^2)

  z <- (kx %*% t(ky)) / (2 * pi * h[1] * h[2] * length(x))
  list(x = gx, y = gy, z = z, h = h)
}

#' Plot 1D Mixture Fit
#'
#' Visualizes a marginal 1D density implied by a \code{geoess_mix} object.
#'
#' @param mix A \code{geoess_mix} mixture object.
#' @param dim_index Integer index of the dimension to plot.
#' @param xlim Optional x-axis limits.
#' @param n Integer grid size.
#' @param components Logical; plot component densities if TRUE.
#' @param add Logical; add to existing plot if TRUE.
#' @param col Line color for the mixture density.
#' @param lwd Line width for the mixture density.
#' @param lty Line type for the mixture density.
#' @param comp_col Line color for component densities.
#' @param comp_lty Line type for component densities.
#' @param main Plot title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param ... Passed to plotting functions.
#'
#' @return Invisible list with grid and densities.
#' @export
geoess_plot_mix_1d <- function(mix, dim_index = 1, xlim = NULL, n = 300,
                               components = TRUE, add = FALSE, col = "black",
                               lwd = 2, lty = 1, comp_col = "gray60",
                               comp_lty = 3, main = NULL, xlab = NULL,
                               ylab = "density", ...) {
  .geoess_check_mix(mix)
  if (!is.numeric(dim_index) || length(dim_index) != 1L || dim_index < 1) {
    stop("dim_index must be a positive integer")
  }
  dim_index <- as.integer(dim_index)
  if (!is.numeric(n) || length(n) != 1L || n <= 10) {
    stop("n must be an integer > 10")
  }

  pars <- .geoess_mix_marginals(mix, dim_index)
  mu_vec <- as.numeric(pars$mu)
  var_vec <- as.numeric(pars$var)
  if (is.null(xlim)) {
    xlim <- .geoess_mix_xlim(mu_vec, var_vec)
  }
  x <- seq(xlim[1], xlim[2], length.out = n)

  dens <- numeric(length(x))
  comp <- matrix(0, nrow = length(x), ncol = length(mu_vec))
  for (k in seq_along(mu_vec)) {
    sd_k <- sqrt(var_vec[k])
    comp[, k] <- pars$w[k] * stats::dnorm(x, mean = mu_vec[k], sd = sd_k)
    dens <- dens + comp[, k]
  }

  if (is.null(xlab)) xlab <- paste0("dim ", dim_index)

  if (!add) {
    plot(x, dens, type = "l", col = col, lwd = lwd, lty = lty,
         xlab = xlab, ylab = ylab, main = main, ...)
  } else {
    lines(x, dens, col = col, lwd = lwd, lty = lty, ...)
  }

  if (components) {
    for (k in seq_len(ncol(comp))) {
      lines(x, comp[, k], col = comp_col, lty = comp_lty, lwd = 1)
    }
  }

  invisible(list(x = x, density = dens, components = comp, dim_index = dim_index))
}

#' Plot 2D Mixture Fit
#'
#' Visualizes a 2D marginal density implied by a \code{geoess_mix} object.
#'
#' @param mix A \code{geoess_mix} mixture object.
#' @param dims Integer vector of length 2 selecting dimensions.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @param n Integer grid size for density evaluation.
#' @param method Character. "contour" or "image".
#' @param add Logical; add to existing plot if TRUE.
#' @param col Colors for image (if method = "image").
#' @param levels Integer number of contour levels.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param ... Passed to plotting functions.
#'
#' @return Invisible list with grid and density matrix.
#' @export
geoess_plot_mix_2d <- function(mix, dims = c(1, 2), xlim = NULL, ylim = NULL,
                               n = 60, method = c("contour", "image"),
                               add = FALSE, col = NULL, levels = 10,
                               xlab = NULL, ylab = NULL, ...) {
  .geoess_check_mix(mix)
  method <- match.arg(method)
  if (!is.numeric(dims) || length(dims) != 2L || any(dims < 1)) {
    stop("dims must be a length-2 integer vector")
  }
  dims <- as.integer(dims)
  if (!is.numeric(n) || length(n) != 1L || n <= 10) {
    stop("n must be an integer > 10")
  }

  pars <- .geoess_mix_marginals(mix, dims)
  if (is.null(xlim)) xlim <- .geoess_mix_xlim(pars$mu[, 1], pars$var[, 1])
  if (is.null(ylim)) ylim <- .geoess_mix_xlim(pars$mu[, 2], pars$var[, 2])

  grid_x <- seq(xlim[1], xlim[2], length.out = n)
  grid_y <- seq(ylim[1], ylim[2], length.out = n)
  grid <- expand.grid(x = grid_x, y = grid_y)

  dens <- numeric(nrow(grid))
  for (k in seq_along(pars$w)) {
    dens <- dens + pars$w[k] * mvtnorm::dmvnorm(
      as.matrix(grid),
      mean = pars$mu[k, ],
      sigma = pars$Sigma[[k]],
      log = FALSE
    )
  }

  z <- matrix(dens, nrow = length(grid_x), ncol = length(grid_y))

  if (is.null(xlab)) xlab <- paste0("dim ", dims[1])
  if (is.null(ylab)) ylab <- paste0("dim ", dims[2])

  if (method == "contour") {
    contour(grid_x, grid_y, z, xlab = xlab, ylab = ylab,
            add = add, nlevels = levels, ...)
  } else {
    if (is.null(col)) col <- gray.colors(20, start = 1, end = 0)
    image(grid_x, grid_y, z, xlab = xlab, ylab = ylab,
          add = add, col = col, ...)
  }

  invisible(list(x = grid_x, y = grid_y, z = z, dims = dims))
}

.geoess_check_mix <- function(mix) {
  if (!inherits(mix, "geoess_mix")) stop("mix must be a geoess_mix object")
  if (is.null(mix$weights) || is.null(mix$mu) || is.null(mix$Sigma)) {
    stop("mix object is missing weights, mu, or Sigma")
  }
  invisible(TRUE)
}

.geoess_mix_marginals <- function(mix, dims) {
  d <- length(mix$mu[[1]])
  if (any(dims > d)) stop("dims exceed mixture dimension")

  K <- length(mix$weights)
  mu <- matrix(NA_real_, nrow = K, ncol = length(dims))
  var <- matrix(NA_real_, nrow = K, ncol = length(dims))
  Sigma <- vector("list", K)

  for (k in seq_len(K)) {
    mu[k, ] <- mix$mu[[k]][dims]
    Sigma_k <- .geoess_as_matrix(mix$Sigma[[k]], "Sigma")
    Sigma_k <- Sigma_k[dims, dims, drop = FALSE]
    Sigma[[k]] <- Sigma_k
    var[k, ] <- diag(Sigma_k)
  }

  list(w = mix$weights, mu = mu, var = var, Sigma = Sigma)
}

.geoess_mix_xlim <- function(mu, var, scale = 4) {
  sd <- sqrt(var)
  rng <- range(c(mu - scale * sd, mu + scale * sd))
  if (!all(is.finite(rng))) stop("Unable to compute xlim from mixture")
  rng
}
