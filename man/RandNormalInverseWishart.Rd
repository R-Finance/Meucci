\name{RandNormalInverseWishart}
\alias{RandNormalInverseWishart}
\title{Generates a multivariate i.i.d. sample of lenght J from the normal-inverse-Wishart distribution, as described in
A. Meucci "Risk and Asset Allocation", Springer, 2005.}
\usage{
  RandNormalInverseWishart(Mu_0, T_0, Sigma_0, nu_0, J)
}
\arguments{
  \item{Mu_0}{: [vector]}

  \item{T_0}{: [scalar]}

  \item{Sigma_0}{: [matrix]}

  \item{nu_0}{: [scalar]}

  \item{J}{: [scalar]}
}
\value{
  Mu : [vector]

  Sigma : [matrix]

  InvSigma : [matrix]
}
\description{
  Generates a multivariate i.i.d. sample of lenght J from
  the normal-inverse-Wishart distribution, as described in
  A. Meucci "Risk and Asset Allocation", Springer, 2005.
}
\note{
  Mu|Sigma ~ N(Mu_0,Sigma/T_0) inv(Sigma) ~
  W(Nu_0,inv(Sigma_0)/Nu_0)
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "RandNormalInverseWishart.m"
}

