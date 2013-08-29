\name{BlackLittermanFormula}
\alias{BlackLittermanFormula}
\title{This function computes the Black-Litterman formula for the moments of the posterior normal, as described in
A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
\usage{
  BlackLittermanFormula(Mu, Sigma, P, v, Omega)
}
\arguments{
  \item{Mu}{: [vector] (N x 1) prior expected values.}

  \item{Sigma}{: [matrix] (N x N) prior covariance matrix.}

  \item{P}{: [matrix] (K x N) pick matrix.}

  \item{v}{: [vector] (K x 1) vector of views.}

  \item{Omega}{: [matrix] (K x K) matrix of confidence.}
}
\value{
  BLMu : [vector] (N x 1) posterior expected values.

  BLSigma : [matrix] (N x N) posterior covariance matrix.
}
\description{
  This function computes the Black-Litterman formula for
  the moments of the posterior normal, as described in A.
  Meucci, "Risk and Asset Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for
  "BlackLittermanFormula.m"
}

