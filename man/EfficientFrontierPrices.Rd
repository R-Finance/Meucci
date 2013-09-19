\name{EfficientFrontierPrices}
\alias{EfficientFrontierPrices}
\title{Computes the mean-variance efficient frontier (on prices) by quadratic programming}
\usage{
  EfficientFrontierPrices(NumPortf, Covariance,
    ExpectedValues, Current_Prices, Budget)
}
\arguments{
  \item{NumPortf}{[scalar] number of portfolio in the
  efficient frontier}

  \item{Covariance}{[matrix] (N x N) covariance matrix}

  \item{ExpectedValues}{[vector] (N x 1) expected returns}

  \item{Current_Prices}{[vector] (N x 1) current prices}

  \item{Budget}{[scalar] budget constraint}
}
\value{
  ExpectedValue [vector] (NumPortf x 1) expected values of
  the portfolios

  Std_Deviation [vector] (NumPortf x 1) standard deviations
  of the portfolios

  Composition [matrix] (NumPortf x N) optimal portfolios
}
\description{
  Compute the mean-variance efficient frontier (on prices)
  by quadratic programming, as described in A. Meucci "Risk
  and Asset Allocation", Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}.

  See Meucci's script for "EfficientFrontierReturns.m".
}

