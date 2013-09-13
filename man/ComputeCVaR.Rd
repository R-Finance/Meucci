\name{ComputeCVaR}
\alias{ComputeCVaR}
\title{Computes the conditional value at risk as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice",
The Risk Magazine, October 2008, p 100-106}
\usage{
  ComputeCVaR(Units, Scenarios, Conf)
}
\arguments{
  \item{Units}{panel of joint factors realizations}

  \item{Scenarios}{vector of probabilities}

  \item{Conf}{Confidence}
}
\value{
  CVaR Conditional Value at Risk
}
\description{
  Computes the conditional value at risk as it appears in
  A. Meucci, "Fully Flexible Views: Theory and Practice",
  The Risk Magazine, October 2008, p 100-106
}
\author{
  Ram Ahluwalia \email{ram@wingedfootcapital.com}
}
\references{
  A. Meucci, "Fully Flexible Views: Theory and Practice"
  \url{http://www.symmys.com/node/158} See Meucci script
  for "ButterflyTrading/ComputeCVaR.m"
}

