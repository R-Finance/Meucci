\name{PlotResults}
\alias{PlotResults}
\title{Plots the results of computing the efficient frontier (Expected returns and frontier), as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine,
October 2008, p 100-106.}
\usage{
  PlotResults(e, s, w, M, Lower = NULL, Upper = NULL)
}
\arguments{
  \item{e}{the NumPortf x 1 matrix of expected returns for
  each portfolio along the efficient frontier}

  \item{s}{the NumPortf x 1 matrix of standard deviation of
  returns for each portfolio along the efficient frontier}

  \item{w}{the NumPortf x N matrix of compositions
  (security weights) for each portfolio along the efficient
  frontier}

  \item{M}{the NumPortf x 1 vector of expected returns for
  each asset}

  \item{Lower}{constraints}

  \item{Upper}{constraints}
}
\description{
  Plots the results of computing the efficient frontier
  (Expected returns and frontier), as it appears in A.
  Meucci, "Fully Flexible Views: Theory and Practice", The
  Risk Magazine, October 2008, p 100-106.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "RankingInformation/PlotResults.m"
}

