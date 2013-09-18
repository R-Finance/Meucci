\name{FitExpectationMaximization}
\alias{FitExpectationMaximization}
\title{Expectation-Maximization (EM) algorithm to recover missing observations in a time series.}
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
  A. Meucci, "Risk and Asset Allocation", Springer, 2005,
  section 4.6.2 "Missing data".
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 177 -
  Expectation-Maximization algorithm for missing data:
  formulas" See Meucci's script for
  "FitExpectationMaximization.m"

  Dempster, A. P. and Laird, M. N. and Rubin, D. B. -
  "Maximum Likelihood from Incomplete Data Via the EM
  Algorithm", Journal of the Royal Statistical Society,
  1977 vol 39 pag. 1-22.

  Bilmes, J. A.- "A Gentle Tutorial of the EM Algorithm and
  its Application to Parameter Estimation for Gaussian
  Mixture and Hidden Markov Models", 1998.
}

