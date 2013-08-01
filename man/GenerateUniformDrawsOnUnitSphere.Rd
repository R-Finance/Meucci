\name{GenerateUniformDrawsOnUnitSphere}
\alias{GenerateUniformDrawsOnUnitSphere}
\title{Generate a uniform sample on the unit hypersphere, as described in  A. Meucci,
 "Risk and Asset Allocation", Springer, 2005.}
\usage{
  GenerateUniformDrawsOnUnitSphere(J, N)
}
\arguments{
  \item{J}{: [scalar] number of draws}

  \item{N}{: [scalar] dimension}
}
\value{
  X : [matrix] (T x N) of draws
}
\description{
  Generate a uniform sample on the unit hypersphere, as
  described in A. Meucci, "Risk and Asset Allocation",
  Springer, 2005.
}
\note{
  Initial script by Xiaoyu Wang - Dec 2006 We decompose
  X=U*R, where U is a uniform distribution on unit sphere
  and
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "GenerateUniformDrawsOnUnitSphere.m"
}

