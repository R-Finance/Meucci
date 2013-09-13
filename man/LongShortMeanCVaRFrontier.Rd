\name{LongShortMeanCVaRFrontier}
\alias{LongShortMeanCVaRFrontier}
\title{Computes the long-short conditional value at risk frontier as it appears in A. Meucci,
"Fully Flexible Views: Theory and Practice", The Risk Magazine, October 2008, p 100-106}
\usage{
  LongShortMeanCVaRFrontier(PnL, Probs, Butterflies,
    Options)
}
\arguments{
  \item{PnL}{Profit and Loss scenarios}

  \item{Probs}{vector of probabilities}

  \item{Butterflies}{list of securities with some analytics
  computed.}

  \item{Options}{list of options}
}
\value{
  Exp vector of expected returns for each asset

  SDev vector of security volatilities along the efficient
  frontier

  CVaR Conditional Value at Risk for each portfolio

  Composition matrix of compositions (security weights) for
  each portfolio along the efficient frontier
}
\description{
  Computes the long-short conditional value at risk
  frontier as it appears in A. Meucci, "Fully Flexible
  Views: Theory and Practice", The Risk Magazine, October
  2008, p 100-106
}
\author{
  Ram Ahluwalia \email{ram@wingedfootcapital.com}, Xavier
  Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "ButterflyTrading/LongShortMeanCVaRFrontier.m"
}

