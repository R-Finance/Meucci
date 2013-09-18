\name{FitMultivariateGarch}
\alias{FitMultivariateGarch}
\title{Estimation of multivariate GARCH models}
\usage{
  FitMultivariateGarch(returns, demean = 1, eps = 0,
    df = 500)
}
\arguments{
  \item{returns}{: [matrix] (T x N) returns so rows must
  correspond to time and columns to assets}

  \item{demean}{: [scalar] specifies whether returns should
  be demeaned (if demean = 1) or not to estimate the model;
  default value is 1.}

  \item{eps}{: [scalar] used in enforcing a_ii + b_ii <= 1
  - eps; the default value is zero}

  \item{df}{: [scalar] degree of freedom for the
  t-distribution; the default value is 500 to make it,
  basically, normal}
}
\value{
  mu : [vector]

  ATMF : [matrix] coefficient matrix A-tilde (in the
  notation of the paper)

  BTMF : [matrix] coefficient matrix B-tilde (in the
  notation of the paper)

  CTMF : [matrix] coefficient matrix C-tilde (in the
  notation of the paper)

  Hhat : [matrix] forecasted conditional covariance matrix
}
\description{
  Estimation of multivariate GARCH models
}
\note{
  Code for MATLAB initially written by Olivier Ledoit and
  Michael Wolf
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 136 -
  Equity market: multivariate GARCH process".

  See Meucci's script for "FitMultivariateGarch.m"
}

