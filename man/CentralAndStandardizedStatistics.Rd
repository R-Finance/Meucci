\name{CentralAndStandardizedStatistics}
\alias{CentralAndStandardizedStatistics}
\title{Compute central and standardized statistics.}
\usage{
  CentralAndStandardizedStatistics(X, N)
}
\arguments{
  \item{X}{: [vector] (J x 1) draws from the distribution}

  \item{N}{: [scalar] highest degree for the central
  moment}
}
\value{
  ga : [vector] (1 x N) standardized statistics up to order
  N

  mu : [vector] (1 x N) central moments up to order N
}
\description{
  Compute central and standardized statistics, as described
  in A. Meucci "Risk and Asset Allocation", Springer, 2005.

  Computes the central moments \deqn{ CM_1^X \equiv
  \mu_{X}\,, \quad CM_n^X \equiv E \{(X - E\{ X
  \})^{n}\}\,, \quad n=2,3,\ldots ,} and from them the
  standarized statistics \deqn{
  \mu_{X},\sigma_{X},sk_{X},ku_{X},\gamma_{X}^{(5)}, \ldots
  ,\gamma_{X}^{(n)} .} where \deqn{\gamma_{X}^{(n)} \equiv
  E \{(X - \mu_{X})^{n}\}/\sigma_{X}^{n},\quad n\geq3 .}
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 97 -
  Projection of skewness, kurtosis, and all standardized
  summary statistics". See Meucci's script for
  "CentralAndStandardizedStatistics.m"

  Kendall, M., Stuart, A., 1969. The Advanced Theory of
  Statistics, Volume, 3rd Edition. Griffin.

  A. Meucci - "Annualization and general projection of
  skweness, kurtosis, and all summary statistics", GARP
  Risk Professional August 2010, 55-56.
  \url{http://symmys.com/node/136}.
}

