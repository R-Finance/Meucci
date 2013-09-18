\name{SimulateJumpDiffusionMerton}
\alias{SimulateJumpDiffusionMerton}
\title{Simulates a Merton jump-diffusion process.}
\usage{
  SimulateJumpDiffusionMerton(m, s, l, a, D, ts, J)
}
\arguments{
  \item{m}{[scalar] deterministic drift of diffusion}

  \item{s}{[scalar] standard deviation of diffusion}

  \item{l}{[scalar] Poisson process arrival rate}

  \item{a}{[scalar] drift of log-jump}

  \item{D}{[scalar] st.dev of log-jump}

  \item{ts}{[vector] time steps}

  \item{J}{[scalar] number of simulations}
}
\value{
  X [matrix] (J x length(ts)) of simulations
}
\description{
  This function simulates a jump diffusion process, as
  described in A. Meucci "Risk and Asset Allocation",
  Springer, 2005.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 132 -
  Simulation of a jump-diffusion process".

  See Meucci's script for "SimulateJumpDiffusionMerton.m"

  Merton, R. C., 1976. "Option pricing when underlying
  stocks are discontinuous". Journal of Financial Economics
  3, 125-144.
}

