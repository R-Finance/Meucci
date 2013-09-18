\name{LognormalMoments2Parameters}
\alias{LognormalMoments2Parameters}
\title{Computes the mean and standard deviation of a lognormal distribution from its parameters.}
\usage{
  LognormalMoments2Parameters(e, v)
}
\arguments{
  \item{e}{[scalar] expected value of the lognormal
  distribution}

  \item{v}{[scalar] variance of the lognormal distribution}
}
\value{
  mu [scalar] expected value of the normal distribution

  sig2 [scalar] variance of the normal distribution
}
\description{
  Computes the mean and standard deviation of a lognormal
  distribution from its parameters, as described in A.
  Meucci, "Risk and Asset Allocation", Springer, 2005.

  \deqn{\sigma^{2} = \ln \left( 1 + \frac{V}{E^{2}} \right)
  , } \deqn{\mu = \ln(E) - \frac{1}{2} \ln \left( 1 +
  \frac{V}{E^{2}} \right) .}
}
\note{
  Inverts the formulas (1.98)-(1.99) in "Risk and Asset
  Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 25 -
  Simulation of a lognormal random variable".

  See Meucci's script for "LognormalMoments2Parameters.m"
}

