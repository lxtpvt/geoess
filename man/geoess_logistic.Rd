% Generated manually to mirror roxygen comments in R/beta_binomial_canonical.R
\name{geoess_logistic}
\alias{geoess_logistic}
\title{Logistic Link Helper}
\usage{
geoess_logistic(eta)
}
\arguments{
\item{eta}{Numeric vector.}
}
\value{
Numeric vector on \eqn{(0, 1)}.
}
\description{
Stable logistic transformation
\eqn{\theta = \operatorname{logit}^{-1}(\eta) = 1 / (1 + e^{-\eta})}.
}
\examples{
geoess_logistic(c(-5, 0, 5))
}
