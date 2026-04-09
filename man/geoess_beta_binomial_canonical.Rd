% Generated manually to mirror roxygen comments in R/beta_binomial_canonical.R
\name{geoess_beta_binomial_canonical}
\alias{geoess_beta_binomial_canonical}
\title{Beta Prior ESS for Binomial/Bernoulli on the Canonical Logit Scale}
\usage{
geoess_beta_binomial_canonical(a, b, n = 1L, return = c("ess", "details"))
}
\arguments{
\item{a}{Numeric shape1 parameter (> 0).}

\item{b}{Numeric shape2 parameter (> 0).}

\item{n}{Optional sample size. Included for context; the returned ESS is
per-trial and therefore does not depend on \code{n}.}

\item{return}{Character. \code{"ess"} returns the scalar ESS. \code{"details"}
returns the ESS together with callable Fisher-information and curvature
functions on the \eqn{\eta}-scale.}
}
\value{
If \code{return = "ess"}, a numeric scalar \eqn{a+b}. If
\code{return = "details"}, a list with \code{ess_local},
\code{ess_global}, \code{parameterization}, \code{I1_eta},
\code{Ipi_eta}, \code{ratio_eta}, and \code{notes}.
}
\description{
Computes the effective sample size (ESS) of a Beta prior for a
Bernoulli/Binomial likelihood on the canonical parameter scale
\eqn{\eta = \operatorname{logit}(\theta)}.
}
\details{
For one Bernoulli observation, the log likelihood in canonical form is
\deqn{\ell(\eta; y) = y\eta - \log(1 + e^\eta),}
with \eqn{\theta = \operatorname{logit}^{-1}(\eta)}. The per-trial Fisher
information on the \eqn{\eta}-scale is
\deqn{I_1(\eta) = \theta(\eta)\{1-\theta(\eta)\}.}

If \eqn{\theta \sim \operatorname{Beta}(a, b)}, then using
\eqn{d\theta/d\eta = \theta(1-\theta)}, the induced prior density on
\eqn{\eta} satisfies
\deqn{\pi(\eta) \propto \theta(\eta)^a \{1-\theta(\eta)\}^b,}
so that
\deqn{\log \pi(\eta) = a \log \theta(\eta) + b \log \{1-\theta(\eta)\} + C.}
Differentiating twice yields
\deqn{I_\pi(\eta) = (a+b)\theta(\eta)\{1-\theta(\eta)\}.}

Therefore the local information ratio is constant:
\deqn{\frac{I_\pi(\eta)}{I_1(\eta)} = a+b,}
so the local ESS and global ESS are both \eqn{a+b} on the canonical scale.
This matches the classical pseudo-count interpretation and the ELIR logic
when expressed on the natural parameter scale.

The argument \code{n} is included to emphasize that the function returns a
\emph{per-trial} ESS. For \eqn{n} Bernoulli trials, the total likelihood
information is \eqn{n I_1(\eta)}, but the ESS relative to one trial remains
\eqn{a+b}.
}
\examples{
geoess_beta_binomial_canonical(3, 7)
# returns 10

det <- geoess_beta_binomial_canonical(3, 7, return = "details")
det$ratio_eta(c(-5, -1, 0, 1, 5))
}
