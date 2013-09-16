\name{MvnRnd}
\alias{MvnRnd}
\title{Generate normal simulations whose sample moments match the population moments}
\usage{
  MvnRnd(M, S, J)
}
\arguments{
  \item{M}{: [vector] (N x 1) expectation}

  \item{S}{: [matrix] (N x N) covariance matrix}

  \item{J}{: [scalar] number of draws (even number)}
}
\value{
  X : [matrix] (J x N) of drawsF_U : [vector] (J x 1) PDF
  values
}
\description{
  Generate normal simulations whose sample moments match
  the population moments, as described in A. Meucci, "Risk
  and Asset Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 64 -
  Simulation of a multivariate normal random variable with
  matching moments".

  See Meucci's script for "MvnRnd.m".
}

