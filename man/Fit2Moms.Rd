\name{Fit2Moms}
\alias{Fit2Moms}
\title{Uses Entropy Pooling to compute a double-decay covariance matrix.}
\usage{
  Fit2Moms(X, m, S)
}
\arguments{
  \item{X}{matrix representing the risk drivers.}

  \item{m}{matrix of zeros, representing the expectation of
  the risk drivers.}

  \item{S}{matrix representing the double-decay estimation
  for the correlation matrix of the risk drivers.}
}
\value{
  p list containing the vector of posterior probabilities
  and information about the optimization performance.
}
\description{
  This function uses Entropy Pooling to compute a
  double-decay covariance matrix, as described in A.
  Meucci, "Personalized Risk Management: Historical
  Scenarios with Fully Flexible Probabilities" GARP Risk
  Professional, Dec 2010, p 47-51
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://www.symmys.com/node/150} See Meucci script
  for "S_MainFullFlexProbs.m"
}

