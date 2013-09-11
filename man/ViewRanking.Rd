\name{ViewRanking}
\alias{ViewRanking}
\title{Computes posterior probabilities to view the rankings, as it appears in A. Meucci,
"Fully Flexible Views: Theory and Practice", The Risk Magazine, October 2008, p 100-106.}
\usage{
  ViewRanking(X, p, Lower, Upper)
}
\arguments{
  \item{X}{a vector containing returns for all the asset
  classes}

  \item{p}{a vector containing the prior probability
  values}

  \item{Lower}{a vector of indexes indicating which column
  is lower than the corresponding column number in Upper}

  \item{Upper}{a vector of indexes indicating which column
  is lower than the corresponding column number in Upper}
}
\description{
  Computes posterior probabilities to view the rankings, as
  it appears in A. Meucci, "Fully Flexible Views: Theory
  and Practice", The Risk Magazine, October 2008, p
  100-106.
}
\author{
  Ram Ahluwalia \email{ram@wingedfootcapital.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "RankingInformation/ViewRanking.m"
}

