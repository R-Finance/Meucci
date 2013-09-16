\name{CentralAndStandardizedStatistics}
\alias{CentralAndStandardizedStatistics}
\title{Compute central and standardized statistics, as described in A. Meucci
"Risk and Asset Allocation", Springer, 2005}
\usage{
  CentralAndStandardizedStatistics(X, N)
}
\arguments{
  \item{X}{: [vector] (J x 1) draws from the distribution}

  \item{N}{: [scalar] highest degree for the central
  moment}
}
\value{
  ga : [vector] (1 x N) standardized statistics up to order
  N

  mu : [vector] (1 x N) central moments up to order N
}
\description{
  Compute central and standardized statistics, as described
  in A. Meucci "Risk and Asset Allocation", Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "CentralAndStandardizedStatistics.m"
}

