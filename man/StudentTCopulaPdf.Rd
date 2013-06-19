\name{StudentTCopulaPdf}
\alias{StudentTCopulaPdf}
\title{Pdf of the copula of the Student t distribution at the generic point u in the unit hypercube,
as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
\usage{
  StudentTCopulaPdf(u, nu, Mu, Sigma)
}
\arguments{
  \item{u}{: [vector] (J x 1) grade}

  \item{nu}{: [numerical] degrees of freedom}

  \item{Mu}{: [vector] (N x 1) mean}

  \item{Sigma}{: [matrix] (N x N) scatter}
}
\value{
  F_U : [vector] (J x 1) PDF values
}
\description{
  Pdf of the copula of the Student t distribution at the
  generic point u in the unit hypercube, as described in A.
  Meucci, "Risk and Asset Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for
  "StudentTCopulaPdf.m"
}

