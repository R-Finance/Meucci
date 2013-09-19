\name{Log2Lin}
\alias{Log2Lin}
\title{Maps moments of log-returns to linear returns .}
\usage{
  Log2Lin(Mu, Sigma)
}
\arguments{
  \item{Mu}{[vector] (N x 1)}

  \item{Sigma}{[matrix] (N x N)}
}
\value{
  M [vector] (N x 1)

  S [matrix] (N x N)
}
\description{
  Map moments of log-returns to linear returns, as
  described in A. Meucci, "Risk and Asset Allocation",
  Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}.

  See Meucci's script for "Log2Lin.m"
}

