.geoess_as_matrix <- function(x, name) {
  if (inherits(x, "Matrix")) {
    if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
    if (any(!is.finite(x))) stop(sprintf("%s contains non-finite entries", name))
    if (nrow(x) == 0 || ncol(x) == 0) {
      stop(sprintf("%s must have positive dimensions", name))
    }
    return(x)
  }

  if (is.null(dim(x))) {
    if (length(x) != 1L || !is.numeric(x)) {
      stop(sprintf("%s must be a numeric matrix or scalar", name))
    }
    x <- matrix(as.numeric(x), nrow = 1L, ncol = 1L)
  }

  if (!is.matrix(x)) stop(sprintf("%s must be a matrix", name))
  if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
  if (any(!is.finite(x))) stop(sprintf("%s contains non-finite entries", name))
  if (nrow(x) == 0 || ncol(x) == 0) {
    stop(sprintf("%s must have positive dimensions", name))
  }
  x
}

.geoess_check_square <- function(x, name) {
  if (nrow(x) == 0 || ncol(x) == 0) {
    stop(sprintf("%s must have positive dimensions", name))
  }
  if (nrow(x) != ncol(x)) {
    stop(sprintf("%s must be square", name))
  }
  invisible(TRUE)
}

.geoess_check_same_dim <- function(a, b, name_a, name_b) {
  if (nrow(a) != nrow(b) || ncol(a) != ncol(b)) {
    stop(sprintf("Dimension mismatch between %s and %s", name_a, name_b))
  }
  invisible(TRUE)
}

.geoess_symmetrize <- function(x) {
  0.5 * (x + t(x))
}

.geoess_diag_scalar <- function(val, d) {
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
    stop("val must be a finite numeric scalar")
  }
  if (!is.numeric(d) || length(d) != 1L || d < 1L) {
    stop("d must be a positive integer")
  }
  d <- as.integer(d)
  m <- matrix(0, nrow = d, ncol = d)
  diag(m) <- as.numeric(val)
  m
}

.geoess_solve <- function(A, b = NULL, jitter = 1e-8) {
  use_matrix <- inherits(A, "Matrix") || (!is.null(b) && inherits(b, "Matrix"))
  solve_fun <- if (use_matrix) Matrix::solve else base::solve

  tryCatch(
    if (is.null(b)) solve_fun(A) else solve_fun(A, b),
    error = function(e) {
      warning("solve failed; adding jitter to diagonal")
      A_j <- if (use_matrix) {
        A + Matrix::Diagonal(n = nrow(A), x = jitter)
      } else {
        A + diag(jitter, nrow(A))
      }
      if (is.null(b)) solve_fun(A_j) else solve_fun(A_j, b)
    }
  )
}

.geoess_try_chol <- function(A, jitter = 1e-8) {
  A_sym <- .geoess_symmetrize(A)
  R <- tryCatch(chol(A_sym), error = function(e) NULL)
  if (!is.null(R)) {
    return(list(R = R, A = A_sym, jittered = FALSE))
  }

  A_j <- A_sym + diag(jitter, nrow(A_sym))
  R <- tryCatch(chol(A_j), error = function(e) NULL)
  list(R = R, A = A_j, jittered = TRUE)
}

.geoess_prec_from_chol <- function(R) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) {
    stop("R must be a square upper-triangular matrix")
  }
  d <- nrow(R)
  A <- backsolve(R, diag(d))
  A %*% t(A)
}

.geoess_trace <- function(A) {
  sum(diag(A))
}

.geoess_psd_adjust <- function(M, psd) {
  M <- .geoess_symmetrize(M)
  if (psd == "none") return(M)
  if (psd == "clip") {
    ev <- eigen(M, symmetric = TRUE)
    ev$values[ev$values < 0] <- 0
    return(ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors))
  }
  as.matrix(Matrix::nearPD(M, ensureSymmetry = TRUE)$mat)
}

.geoess_kl_normal_var <- function(v1, v2) {
  0.5 * (log(v2 / v1) + (v1 / v2) - 1)
}

.geoess_kl_discrete <- function(p, q, eps = 1e-12) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  if (length(p) != length(q)) stop("p and q must be the same length")
  if (any(p < 0) || any(q < 0)) stop("p and q must be nonnegative")

  p <- p / sum(p)
  q <- pmax(q, eps)
  q <- q / sum(q)

  idx <- p > 0
  sum(p[idx] * log(p[idx] / q[idx]))
}
