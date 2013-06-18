\name{LognormalParam2Statistics}
\alias{LognormalParam2Statistics}
\title{Compute expectation, Cov, standard deviation and Corr for a lognormal distribution, as described in
A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  LognormalParam2Statistics(Mu, Sigma)
}
\arguments{
  \item{Mu}{: [vector] (N x 1) location parameter}

  \item{Sigma}{: [matrix] (N x N) scale parameter}
}
\value{
  Exp : [vector] (N x 1) expectation

  Cov : [matrix] (N x N) covariance

  Std : [vector] (N x 1) standard deviation

  Corr : [matrix] (N x N) correlation
}
\description{
  Compute expectation, Cov, standard deviation and Corr for
  a lognormal distribution, as described in A. Meucci "Risk
  and Asset Allocation", Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for
  "LognormalParam2Statistics.m"
}

