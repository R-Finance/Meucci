\name{ConvertChangeInYield2Price}
\alias{ConvertChangeInYield2Price}
\title{Convert change in yield-to-maturity to price for fixed-income securities, as described in
A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  ConvertChangeInYield2Price(Exp_DY, Cov_DY, Times2Mat,
    CurrentPrices)
}
\arguments{
  \item{Exp_DY}{: [vector] (N x 1) expected value of change
  in yield to maturity}

  \item{Cov_DY}{: [matrix] (N x N) covariance of change in
  yield to maturity}

  \item{Times2Mat}{: [scalar] time to maturity}

  \item{CurrentPrices}{: [vector] (N x 1) current prices}
}
\value{
  Exp_Prices : [vector] (N x 1) expected prices

  Cov_Prices : [matrix] (N x N) covariance of prices
}
\description{
  Convert change in yield-to-maturity to price for
  fixed-income securities, as described in A. Meucci "Risk
  and Asset Allocation", Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See (6.77)-(6.79) in "Risk and Asset
  Allocation"-Springer (2005), by A. Meucci See Meucci's
  script for "ConvertChangeInYield2Price.m"
}

