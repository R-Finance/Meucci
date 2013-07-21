\name{PlotMarginalsNormalInverseWishart}
\alias{PlotMarginalsNormalInverseWishart}
\title{Plot the marginals of the normal-inverse-Whishart model.
Described in A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  PlotMarginalsNormalInverseWishart(Mu_Simul,
    InvSigma_Simul, Mu_0, T_0, Sigma_0, Nu_0, Legend)
}
\arguments{
  \item{Mu_Simul}{: []}

  \item{InvSigma_Simul}{: []}

  \item{Mu_0}{: []}

  \item{T_0}{: []}

  \item{Sigma_0}{: []}

  \item{Nu_0}{: []}

  \item{Legend}{: []}
}
\description{
  Plot the marginals of the normal-inverse-Whishart model.
  Described in A. Meucci "Risk and Asset Allocation",
  Springer, 2005
}
\note{
  Numerically and analytically the marginal pdf of - the
  first entry of the random vector Mu - the (1,1)-entry of
  the random matrix inv(Sigma) when Mu and Sigma are
  jointly normal-inverse-Wishart: Mu ~ St(Mu_0,Sigma/T_0)
  inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "QuantileMixture.m"
}

