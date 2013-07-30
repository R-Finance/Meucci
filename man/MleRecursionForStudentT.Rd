\name{MleRecursionForStudentT}
\alias{MleRecursionForStudentT}
\title{Compute recursively the ML estimators of location and scatter of a multivariate Student t distribution with
given degrees of freedom, as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
\usage{
  MleRecursionForStudentT(x, Nu, Tolerance = 10^(-10))
}
\arguments{
  \item{x}{: [matrix] (T x N) observations}

  \item{Nu}{: [scalar] degrees of freedom parameter}

  \item{Tolerance}{: [scalar] tolerance parameter. Default:
  10^(-10)}
}
\value{
  Mu : [vector] (N x 1) mean

  Sigma : [matrix] (N x N) covariance
}
\description{
  Compute recursively the ML estimators of location and
  scatter of a multivariate Student t distribution with
  given degrees of freedom, as described in A. Meucci,
  "Risk and Asset Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for
  "MleRecursionForStudentT.m"
}

