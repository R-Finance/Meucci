\name{PlotFrontier}
\alias{PlotFrontier}
\title{Plots the efficient frontier, as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine,
 October 2008, p 100-106.}
\usage{
  PlotFrontier(e, s, w)
}
\arguments{
  \item{e}{the NumPortf x 1 matrix of expected returns for
  each portfolio along the efficient frontier}

  \item{s}{the NumPortf x 1 matrix of standard deviation of
  returns for each portfolio along the efficient frontier}

  \item{w}{the NumPortf x N matrix of compositions
  (security weights) for each portfolio along the efficient
  frontier}
}
\description{
  Plots the efficient frontier, as it appears in A. Meucci,
  "Fully Flexible Views: Theory and Practice", The Risk
  Magazine, October 2008, p 100-106.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "RankingInformation/PlotFrontier.m"
}

