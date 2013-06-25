\name{ProjectionStudentT}
\alias{ProjectionStudentT}
\title{Perform the horizon projection of a Student t invariant, as described in
A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  ProjectionStudentT(nu, m, s, T)
}
\arguments{
  \item{nu}{: [scalar] degree of freedom}

  \item{s}{: [scalar] scatter parameter}

  \item{m}{: [scalar] location parameter}

  \item{T}{: [scalar] multiple of the estimation period to
  the invesment horizon}
}
\value{
  x_Hor : [scalar]

  f_Hor : [scalar]

  F_Hor : [scalar]
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
  \url{http://} See Meucci's script for
  "ProjectionStudentT.m"
}

