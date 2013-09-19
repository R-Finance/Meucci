\name{ConvertCompoundedReturns2Price}
\alias{ConvertCompoundedReturns2Price}
\title{Convert compounded returns to prices for equity-like securities.}
\usage{
  ConvertCompoundedReturns2Price(Exp_Comp_Rets,
    Cov_Comp_Rets, Starting_Prices)
}
\arguments{
  \item{Exp_Comp_Rets}{[vector] (N x 1) expected values of
  compounded returns}

  \item{Cov_Comp_Rets}{[matrix] (N x N) covariance matrix
  of compounded returns}

  \item{Starting_Prices}{[vector] (N x 1)}
}
\value{
  Exp_Prices [vector] (N x 1) expected values of prices

  Cov_Prices [matrix] (N x N) covariance matrix of prices
}
\description{
  Convert compounded returns to prices for equity-like
  securities, as described in A. Meucci "Risk and Asset
  Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "ConvertCompoundedReturns2Price.m".

  A. Meucci - "Risk and Asset Allocation"-Springer (2005).
  See (6.77)-(6.79).
}

