\name{NormalCopulaPdf}
\alias{NormalCopulaPdf}
\title{Computes the pdf of the copula of the normal distribution at the generic point u in the unit hypercube}
\usage{
  NormalCopulaPdf(u, Mu, Sigma)
}
\arguments{
  \item{u}{[vector] (J x 1) grade}

  \item{Mu}{[vector] (N x 1) mean}

  \item{Sigma}{[matrix] (N x N) covariance}
}
\value{
  F_U [vector] (J x 1) PDF values
}
\description{
  Computes the pdf of the copula of the normal distribution
  at the generic point u in the unit hypercube, as
  described in A. Meucci, "Risk and Asset Allocation",
  Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 33 - Pdf
  of the normal copula".

  See Meucci's script for "NormalCopulaPdf.m"
}

