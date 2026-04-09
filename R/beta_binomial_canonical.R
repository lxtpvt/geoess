#' Logistic Link Helper
#'
#' Stable logistic transformation
#' \eqn{\theta = \operatorname{logit}^{-1}(\eta) = 1 / (1 + e^{-\eta})}.
#'
#' @param eta Numeric vector.
#'
#' @return Numeric vector on \eqn{(0, 1)}.
#'
#' @examples
#' geoess_logistic(c(-5, 0, 5))
#' @export
geoess_logistic <- function(eta) {
  eta <- as.numeric(eta)
  if (any(!is.finite(eta))) stop("eta must contain only finite values")

  out <- numeric(length(eta))
  pos <- eta >= 0
  out[pos] <- 1 / (1 + exp(-eta[pos]))
  exp_eta <- exp(eta[!pos])
  out[!pos] <- exp_eta / (1 + exp_eta)
  out
}

.geoess_check_beta_shapes <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L) {
    stop("a and b must be numeric scalars")
  }
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("a and b must be positive and finite")
  }
}

.geoess_check_n_trials <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("n must be a positive scalar")
  }
  if (abs(n - round(n)) > sqrt(.Machine$double.eps)) {
    stop("n must be an integer-valued positive scalar")
  }
  as.integer(round(n))
}

#' Beta Prior ESS for Binomial/Bernoulli on the Canonical Logit Scale
#'
#' Computes the effective sample size (ESS) of a Beta prior for a
#' Bernoulli/Binomial likelihood on the canonical parameter scale
#' \eqn{\eta = \operatorname{logit}(\theta)}.
#'
#' For one Bernoulli observation, the log likelihood in canonical form is
#' \deqn{\ell(\eta; y) = y\eta - \log(1 + e^\eta),}
#' with \eqn{\theta = \operatorname{logit}^{-1}(\eta)}. The per-trial Fisher
#' information on the \eqn{\eta}-scale is
#' \deqn{I_1(\eta) = \theta(\eta)\{1-\theta(\eta)\}.}
#'
#' If \eqn{\theta \sim \operatorname{Beta}(a, b)}, then using
#' \eqn{d\theta/d\eta = \theta(1-\theta)}, the induced prior density on
#' \eqn{\eta} satisfies
#' \deqn{\pi(\eta) \propto \theta(\eta)^a \{1-\theta(\eta)\}^b,}
#' so that
#' \deqn{\log \pi(\eta) = a \log \theta(\eta) + b \log \{1-\theta(\eta)\} + C.}
#' Differentiating twice yields the prior curvature
#' \deqn{I_\pi(\eta) = -\frac{d^2}{d\eta^2}\log \pi(\eta)
#' = (a+b)\theta(\eta)\{1-\theta(\eta)\}.}
#'
#' Therefore the local information ratio is constant:
#' \deqn{\frac{I_\pi(\eta)}{I_1(\eta)} = a+b,}
#' so the local ESS and global ESS are both \eqn{a+b} on the canonical scale.
#' This matches the classical pseudo-count interpretation and the ELIR logic
#' when expressed on the natural parameter scale.
#'
#' The argument \code{n} is included to remind users that the function returns
#' a \emph{per-trial} ESS. For \eqn{n} Bernoulli trials, the total likelihood
#' information is \eqn{n I_1(\eta)}, but the ESS relative to one trial remains
#' \eqn{a+b}.
#'
#' @param a Numeric shape1 parameter (> 0).
#' @param b Numeric shape2 parameter (> 0).
#' @param n Optional sample size. Included for context; the returned ESS is
#'   per-trial and therefore does not depend on \code{n}.
#' @param return Character. \code{"ess"} returns the scalar ESS. \code{"details"}
#'   returns the ESS together with callable Fisher-information and curvature
#'   functions on the \eqn{\eta}-scale.
#'
#' @return If \code{return = "ess"}, a numeric scalar \eqn{a+b}. If
#'   \code{return = "details"}, a list with:
#' \itemize{
#'   \item \code{ess_local}: scalar \eqn{a+b}.
#'   \item \code{ess_global}: scalar \eqn{a+b}.
#'   \item \code{parameterization}: \code{"canonical_logit"}.
#'   \item \code{I1_eta}: function computing unit Fisher information
#'   \eqn{\theta(1-\theta)}.
#'   \item \code{Ipi_eta}: function computing prior curvature
#'   \eqn{(a+b)\theta(1-\theta)}.
#'   \item \code{ratio_eta}: function returning the constant ratio \eqn{a+b}.
#'   \item \code{notes}: explanation that the ESS is per-trial on the canonical
#'   scale and does not depend on \code{n}.
#' }
#'
#' @examples
#' geoess_beta_binomial_canonical(3, 7)
#' # returns 10
#'
#' det <- geoess_beta_binomial_canonical(3, 7, return = "details")
#' det$ratio_eta(c(-5, -1, 0, 1, 5))
#' @export
geoess_beta_binomial_canonical <- function(a, b, n = 1L, return = c("ess", "details")) {
  .geoess_check_beta_shapes(a, b)
  n <- .geoess_check_n_trials(n)
  return <- match.arg(return)

  ess_val <- as.numeric(a + b)

  I1_eta <- function(eta) {
    theta <- geoess_logistic(eta)
    theta * (1 - theta)
  }

  Ipi_eta <- function(eta) {
    ess_val * I1_eta(eta)
  }

  ratio_eta <- function(eta) {
    eta <- as.numeric(eta)
    if (any(!is.finite(eta))) stop("eta must contain only finite values")
    rep(ess_val, length(eta))
  }

  if (return == "ess") {
    return(ess_val)
  }

  list(
    ess_local = ess_val,
    ess_global = ess_val,
    parameterization = "canonical_logit",
    n = n,
    I1_eta = I1_eta,
    Ipi_eta = Ipi_eta,
    ratio_eta = ratio_eta,
    notes = paste0(
      "Per-trial ESS on the canonical logit scale for a Beta(", a, ", ", b,
      ") prior under a Bernoulli/Binomial likelihood. ",
      "The unit Fisher information is theta(1-theta), the prior curvature is ",
      "(a+b)theta(1-theta), and the local/global ESS is therefore a+b = ",
      format(ess_val, trim = TRUE), ". For n = ", n,
      " trials, total likelihood information is n * I1_eta(eta), ",
      "but the ESS relative to one trial remains unchanged."
    )
  )
}
