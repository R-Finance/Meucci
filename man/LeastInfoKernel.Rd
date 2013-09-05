\name{LeastInfoKernel}
\alias{LeastInfoKernel}
\title{Computes least information kernel smoothing}
\usage{
  LeastInfoKernel(Y, y, h2)
}
\arguments{
  \item{Y}{Matrix representing the macroeconomic indicator}

  \item{y}{scalar reprenting the target to which Y is
  expected to be close in the Generalized Empirical
  Distribution}

  \item{h2}{N X N matrix}
}
\value{
  p list containing the vector of posterior probabilities
  and information about the optimization performance.
}
\description{
  This script uses Entropy Pooling to compute least
  information kernel smoothing, as described in A. Meucci,
  "Personalized Risk Management: Historical Scenarios with
  Fully Flexible Probabilities" GARP Risk Professional, Dec
  2010, p 47-51
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://www.symmys.com/node/150} See Meucci script
  for "LeastInfoKernel.m"
}

