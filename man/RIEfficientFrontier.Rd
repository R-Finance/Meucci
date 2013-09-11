\name{RIEfficientFrontier}
\alias{RIEfficientFrontier}
\title{Generates an efficient frontier based on Meucci's Ranking Information version and returns a A list with
NumPortf efficient portfolios whos returns are equally spaced along the whole range of the efficient frontier,
as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine, October 2008,
p 100-106.}
\usage{
  RIEfficientFrontier(X, p, Options)
}
\arguments{
  \item{X}{a matrix with the joint-scenario probabilities
  by asset (rows are joint-scenarios, columns are assets)}

  \item{p}{a vector of probabilities associated with each
  scenario in matrix X}

  \item{Options}{a list of options....TBD}
}
\value{
  Exps the NumPortf x 1 vector of expected returns for each
  asset

  Covs the NumPortf x N vector of security volatilities
  along the efficient frontier

  w the NumPortf x N matrix of compositions (security
  weights) for each portfolio along the efficient frontier

  e the NumPortf x 1 matrix of expected returns for each
  portfolio along the efficient frontier

  s the NumPortf x 1 matrix of standard deviation of
  returns for each portfolio along the efficient frontier
}
\description{
  Most recent version of article and MATLAB code available
  at http://www.symmys.com/node/158
}
\author{
  Ram Ahluwalia \email{ram@wingedfootcapital.com} and
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "RankingInformation/EfficientFrontier.m"
}

