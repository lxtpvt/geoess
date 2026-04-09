# Internal helpers for log-scale mixture curvature.

.geoess_Ipi_from_logmix <- function(mix, theta, psd = c("none", "clip", "nearPD"),
                                    eps = 1e-8) {
  psd <- match.arg(psd)
  theta <- as.numeric(theta)
  if (length(theta) != 1L || !is.finite(theta)) return(NULL)
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) return(NULL)

  theta <- max(theta, eps)
  z <- log(theta)
  g <- tryCatch(geoess_mix_gradlog(mix, z), error = function(e) NULL)
  H <- tryCatch(geoess_mix_hesslog(mix, z, symmetrize = TRUE), error = function(e) NULL)
  if (is.null(g) || is.null(H)) return(NULL)
  g <- as.numeric(g)[1]
  H <- as.numeric(H)[1]

  # Transformation: log theta -> theta
  Ipi_val <- -(H - g + 1) / (theta^2)
  if (!is.finite(Ipi_val)) return(NULL)
  if (psd != "none" && Ipi_val < 0) Ipi_val <- 0
  matrix(Ipi_val, nrow = 1L, ncol = 1L)
}

.geoess_samples_mix_log <- function(draws, I1, theta_ref,
                                    psd = c("none", "clip", "nearPD"),
                                    K_grid = 1:2, criterion = "BIC",
                                    eps = 1e-8) {
  psd <- match.arg(psd)
  x <- as.numeric(draws)
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 5L) return(c(local = NA_real_, prior = NA_real_))

  z <- log(pmax(x, eps))
  mix <- tryCatch(
    geoess_automixfit(matrix(z, ncol = 1L), family = "norm",
                      K_grid = K_grid, criterion = criterion),
    error = function(e) NULL
  )
  if (is.null(mix)) return(c(local = NA_real_, prior = NA_real_))

  Ipi_fun <- function(th) .geoess_Ipi_from_logmix(mix, th, psd = psd, eps = eps)
  Ipi_ref <- Ipi_fun(theta_ref)
  if (is.null(Ipi_ref)) return(c(local = NA_real_, prior = NA_real_))
  I1_ref <- if (is.function(I1)) I1(theta_ref) else I1
  ess_local <- tryCatch(geoess_trace(I1_ref, Ipi_ref), error = function(e) NA_real_)
  ess_prior <- tryCatch(
    geoess_global(Q = matrix(x, ncol = 1L), I1 = I1, Ipi = Ipi_fun,
                  method = "trace", psd = psd)$estimate,
    error = function(e) NA_real_
  )
  c(local = ess_local, prior = ess_prior)
}
