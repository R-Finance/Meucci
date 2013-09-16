\name{BlackScholesCallPrice}
\alias{BlackScholesCallPrice}
\alias{BlackScholesCallPutPrice}
\alias{BlackScholesPutPrice}
\title{Compute the Black-Scholes price of a European call or put option
 as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.}
\usage{
  BlackScholesCallPrice(spot, K, r, vol, T)

  BlackScholesPutPrice(spot, K, r, vol, T)

  BlackScholesCallPutPrice(spot, K, r, vol, T)
}
\arguments{
  \item{spot}{: [scalar] spot price of underlying}

  \item{K}{: [scalar] strike of the call optioon}

  \item{r}{: [scalar] risk free rate as a fraction}

  \item{vol}{: [scalar] volatility of the underlying as a
  fraction}

  \item{T}{: [scalar] time to maturity in years}
}
\value{
  c : [scalar] price of European call(s)

  p : [scalar] price of European put(s)

  delta : [scalar] delta of the call(s) or put(s)

  cash : [scalar] cash held in a replicating portfolio
}
\description{
  Compute the Black-Scholes price of a European call or put
  option as described in A. Meucci, "Risk and Asset
  Allocation", Springer, 2005.
}
\note{
  Code is vectorized, so the inputs can be vectors or
  matrices (but sizes must match)
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}. See
  Meucci's script for "BlackScholesCallPrice.m"
}

