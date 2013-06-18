\name{LognormalMoments2Parameters}
\alias{LognormalMoments2Parameters}
\title{Compute the mean and standard deviation of a lognormal distribution from its parameters, as described in
A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 1.}
\usage{
  LognormalMoments2Parameters(e, v)
}
\arguments{
  \item{e}{: [scalar] expected value of the lognormal
  distribution}

  \item{v}{: [scalar] variance of the lognormal
  distribution}
}
\value{
  mu : [scalar] expected value of the normal distribution

  sig2 : [scalar] variance of the normal distribution
}
\description{
  Compute the mean and standard deviation of a lognormal
  distribution from its parameters, as described in A.
  Meucci, "Risk and Asset Allocation", Springer, 2005,
  Chapter 1.
}
\note{
  Inverts the formulas (1.98)-(1.99) in Risk and Asset
  Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for
  "LognormalMoments2Parameters"
}

