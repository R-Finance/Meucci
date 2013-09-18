\name{ProjectionStudentT}
\alias{ProjectionStudentT}
\title{Perform the horizon projection of a Student t invariant}
\usage{
  ProjectionStudentT(nu, m, s, T)
}
\arguments{
  \item{nu}{[scalar] degree of freedom}

  \item{s}{[scalar] scatter parameter}

  \item{m}{[scalar] location parameter}

  \item{T}{[scalar] multiple of the estimation period to
  the invesment horizon}
}
\value{
  x_Hor [scalar] probabilities at horizon

  f_Hor [scalar] horizon discretized pdf (non-standarized)

  F_Hor [scalar] horizon discretized cdf (non-standarized)
}
\description{
  Perform the horizon projection of a Student t invariant,
  as described in A. Meucci "Risk and Asset Allocation",
  Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 141 -
  Fixed-income market: projection of Student t invariants".

  See Meucci's script for "ProjectionStudentT.m"
}

