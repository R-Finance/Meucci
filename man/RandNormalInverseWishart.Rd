\name{RandNormalInverseWishart}
\alias{RandNormalInverseWishart}
\title{Generates a multivariate i.i.d. sample of lenght J from the normal-inverse-Wishart distribution.}
\usage{
  RandNormalInverseWishart(Mu_0, T_0, Sigma_0, nu_0, J)
}
\arguments{
  \item{Mu_0}{[vector] location parameter.}

  \item{T_0}{[scalar] number of observations.}

  \item{Sigma_0}{[matrix] scatter parameter.}

  \item{nu_0}{[scalar] degrees of freedom.}

  \item{J}{[scalar] number of simulations to compute.}
}
\value{
  Mu [vector] location parameter from the
  normal-inverse-Wishart distribution.

  Sigma [matrix] dispersion parameter from the
  normal-inverse-Wishart distribution.

  InvSigma [matrix] inverse of the dispersion parameter
  from the normal-inverse-Wishart distribution.
}
\description{
  Generates a multivariate i.i.d. sample of lenght J from
  the normal-inverse-Wishart distribution, as described in
  A. Meucci "Risk and Asset Allocation", Springer, 2005.
}
\note{
  \deqn{\mu\| \Sigma \sim N(\mu_{0}, \frac{\Sigma}{T_{0}})
  }{Mu|Sigma ~ N(Mu_0,Sigma/T_0)} \deqn{\Sigma^{-1} \sim
  W(\nu_{0},\frac{\Sigma_{0}^{-1}}{\nu_{0}})}{inv(Sigma) ~
  W(Nu_0,inv(Sigma_0)/Nu_0)}
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "RandNormalInverseWishart.m"
}

