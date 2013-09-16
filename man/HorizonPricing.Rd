\name{HorizonPricing}
\alias{HorizonPricing}
\title{Compute the pricing in the horizon, as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice",
 The Risk Magazine, October 2008, p 100-106.}
\usage{
  HorizonPricing(Butterflies, X)
}
\arguments{
  \item{Butterflies}{List of securities with some analytics
  computed.}

  \item{X}{Panel of joint factors realizations}
}
\value{
  PnL Matrix of profit and loss scenarios
}
\description{
  Compute the pricing in the horizon, as it appears in A.
  Meucci, "Fully Flexible Views: Theory and Practice", The
  Risk Magazine, October 2008, p 100-106.
}
\author{
  Ram Ahluwalia \email{ram@wingedfootcapital.com} and
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "ButterflyTrading/HorizonPricing.m"
}

