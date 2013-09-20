\name{pHistPriorPosterior}
\alias{pHistPriorPosterior}
\title{Plot prior and posterior distributions.}
\usage{
  pHistPriorPosterior(X, p, p_)
}
\arguments{
  \item{X}{: [matrix] (J x N) simulations}

  \item{p}{: [vector] (J x 1) prior probabilities}

  \item{p_}{: [vector] (J x 1) posterior probabilities}
}
\description{
  Plot prior and posterior distributions, as described in
  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}.

  See Meucci's script for "pHistPriorPosterior.m"
}

