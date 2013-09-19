\name{EfficientFrontierReturnsBenchmark}
\alias{EfficientFrontierReturnsBenchmark}
\title{Computes the mean-variance efficient frontier (on returns) by quadratic programming.}
\usage{
  EfficientFrontierReturnsBenchmark(NumPortf, Covariance,
    ExpectedValues, Benchmark, Constraints = NULL)
}
\arguments{
  \item{NumPortf}{[scalar] number of portfolio in the
  efficient frontier}

  \item{Covariance}{[matrix] (N x N) covariance matrix}

  \item{ExpectedValues}{[vector] (N x 1) expected returns}

  \item{Benchmark}{[vector] (N x 1) of benchmark weights}

  \item{Constraints}{[struct] set of constraints. Default:
  weights sum to one, and no-short positions}
}
\value{
  ExpectedValue [vector] (NumPortf x 1) expected values of
  the portfolios

  Volatility [vector] (NumPortf x 1) standard deviations of
  the portfolios

  Composition [matrix] (NumPortf x N) optimal portfolios
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
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}.

  See Meucci's script for
  "EfficientFrontierReturnsBenchmark.m"
}

