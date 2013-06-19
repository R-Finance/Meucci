\name{LognormalCopulaPdf}
\alias{LognormalCopulaPdf}
\title{Computes the pdf of the copula of the lognormal distribution at the generic point u in the unit hypercube,
as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
\usage{
  LognormalCopulaPdf(u, Mu, Sigma)
}
\arguments{
  \item{u}{: [vector] (J x 1) grades}

  \item{Mu}{: [vector] (N x 1) location parameter}

  \item{Sigma}{: [matrix] (N x N) scatter parameter}
}
\value{
  F_U : [vector] (J x 1) PDF values
}
\description{
  Computes the pdf of the copula of the lognormal
  distribution at the generic point u in the unit
  hypercube, as described in A. Meucci, "Risk and Asset
  Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for
  "LognormalCopulaPdf.m"
}

