\name{FitExpectationMaximization}
\alias{FitExpectationMaximization}
\title{Expectation-Maximization (EM) algorithm to recover missing observations in a time series ,
as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
\usage{
  FitExpectationMaximization(X)
}
\arguments{
  \item{X}{: [matrix] (T x N) of data}
}
\value{
  E_EM : [vector] (N x 1) expectation

  S_EM : [matrix] (N x N) covariance matrix

  Y : [matrix] (T x N) updated data

  CountLoop : [scalar] number of iterations of the
  algorithm
}
\description{
  Expectation-Maximization (EM) algorithm to recover
  missing observations in a time series , as described in
  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "FitExpectationMaximization.m"
}

