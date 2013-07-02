\name{SimulateJumpDiffusionMerton}
\alias{SimulateJumpDiffusionMerton}
\title{This function simulates a jump diffusion process, as described in A. Meucci "Risk and Asset Allocation",
Springer, 2005}
\usage{
  SimulateJumpDiffusionMerton(m, s, l, a, D, ts, J)
}
\arguments{
  \item{m}{: [scalar] deterministic drift of diffusion}

  \item{s}{: [scalar] standard deviation of diffusion}

  \item{l}{: [scalar] Poisson process arrival rate}

  \item{a}{: [scalar] drift of log-jump}

  \item{D}{: [scalar] st.dev of log-jump}

  \item{ts}{: [vector] time steps}

  \item{J}{: [scalar] number of simulations}
}
\value{
  X : [matrix] (J x length(ts)) of simulations
}
\description{
  This function simulates a jump diffusion process, as
  described in A. Meucci "Risk and Asset Allocation",
  Springer, 2005
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "SimulateJumpDiffusionMerton.m"
}

