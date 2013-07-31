\name{MvnRnd}
\alias{MvnRnd}
\title{Generate normal simulations whose sample moments match the population moments,
as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
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
  Xavier Valls \email{flamejat@gmail.com} and Ram Ahluwalia
  \email{rahluwalia@gmail.com}
}
\references{
  \url{http://symmys.com/node/170},
  \url{http://www.symmys.com/node/162}{A. Meucci -
  "Simulations with Exact Means and Covariances", Risk,
  July 2009} See Meucci's script for "MvnRnd.m"
}

