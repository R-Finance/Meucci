\name{QuantileMixture}
\alias{QuantileMixture}
\title{Computes the quantile of a mixture distirbution by linear interpolation/extrapolation of the cdf.the confidence
level p can be vector. If this vector is uniformly distributed on [0,1] the sample Q is distributed as the mixture.
Described in A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  QuantileMixture(p, a, m_Y, s_Y, m_Z, s_Z)
}
\arguments{
  \item{p}{: [scalar] in [0,1], probability}

  \item{a}{: [scalar] in (0,1), mixing probability}

  \item{m_Y}{: [scalar] mean of normal component}

  \item{s_Y}{: [scalar] standard deviation of normal
  component}

  \item{m_Z}{: [scalar] first parameters of the log-normal
  component}

  \item{s_Z}{: [scalar] second parameter of the log-normal
  component}
}
\value{
  Q : [scalar] quantile
}
\description{
  Computes the quantile of a mixture distirbution by linear
  interpolation/extrapolation of the cdf.the confidence
  level p can be vector. If this vector is uniformly
  distributed on [0,1] the sample Q is distributed as the
  mixture. Described in A. Meucci "Risk and Asset
  Allocation", Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "QuantileMixture.m"
}

