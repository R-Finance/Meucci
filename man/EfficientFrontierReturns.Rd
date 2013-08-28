\name{EfficientFrontierReturns}
\alias{EfficientFrontierReturns}
\title{Compute the mean-variance efficient frontier (on returns) by quadratic programming, as described in
A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  EfficientFrontierReturns(NumPortf, Covariance,
    ExpectedValues, Constraints = NULL)
}
\arguments{
  \item{NumPortf}{: [scalar] number of portfolio in the
  efficient frontier}

  \item{Covariance}{: [matrix] (N x N) covariance matrix}

  \item{ExpectedValues}{: [vector] (N x 1) expected
  returns}

  \item{Constraints}{: [struct] set of constraints.
  Default: weights sum to one, and no-short positions}
}
\value{
  ExpectedValue : [vector] (NumPortf x 1) expected values
  of the portfolios

  Volatility : [vector] (NumPortf x 1) standard deviations
  of the portfolios

  Composition : [matrix] (NumPortf x N) optimal portfolios
}
\description{
  Compute the mean-variance efficient frontier (on returns)
  by quadratic programming, as described in A. Meucci "Risk
  and Asset Allocation", Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "EfficientFrontierReturns.m"
}

