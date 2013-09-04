\name{FitOrnsteinUhlenbeck}
\alias{FitOrnsteinUhlenbeck}
\title{Fit a multivariate OU process at estimation step tau, as described in  A. Meucci
"Risk and Asset Allocation", Springer, 2005}
\usage{
  FitOrnsteinUhlenbeck(Y, tau)
}
\arguments{
  \item{Y}{: [matrix] (T x N)}

  \item{tau}{: [scalar] time step}
}
\value{
  Mu : [vector] long-term means

  Th : [matrix] whose eigenvalues have positive real part /
  mean reversion speed

  Sig : [matrix] Sig = S * S', covariance matrix of
  Brownian motions
}
\description{
  Fit a multivariate OU process at estimation step tau, as
  described in A. Meucci "Risk and Asset Allocation",
  Springer, 2005
}
\note{
  o dY_t = -Th * (Y_t - Mu) * dt + S * dB_t where o dB_t:
  vector of Brownian motions
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "FitOrnsteinUhlenbeck.m"
}

