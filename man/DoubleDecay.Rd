\name{DoubleDecay}
\alias{DoubleDecay}
\title{Computes a double-decay covariance matrix.}
\usage{
  DoubleDecay(X, lmd_c, lmd_s)
}
\arguments{
  \item{X}{matrix representing the risk drivers.}

  \item{lmd_c}{numeric representing the low decay (long
  half-life) for the correlations.}

  \item{lmd_s}{numeric representing the high decay (short
  half-life) for the volatilities.}
}
\value{
  m matrix of zeros, representing the expectation of the
  risk drivers.

  S matrix representing the double-decay estimation for the
  correlation matrix of the risk drivers.
}
\description{
  This function computes a double-decay covariance matrix
  for the risk drivers provided, as described in A. Meucci,
  "Personalized Risk Management: Historical Scenarios with
  Fully Flexible Probabilities" GARP Risk Professional, Dec
  2010, p 47-51
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://www.symmys.com/node/150} See Meucci script
  for "DoubleDecay.m"
}

