\name{MleRecursionForStudentT}
\alias{MleRecursionForStudentT}
\title{Compute recursively the ML estimators of location and scatter of a multivariate Student t distribution with
given degrees of freedom.}
\usage{
  MleRecursionForStudentT(x, Nu, Tolerance = 10^(-10))
}
\arguments{
  \item{x}{: [matrix] (T x N) observations}

  \item{Nu}{: [scalar] degrees of freedom parameter}

  \item{Tolerance}{: [scalar] tolerance parameter. Default:
  10^(-10)}
}
\value{
  Mu : [vector] (N x 1) mean

  Sigma : [matrix] (N x N) covariance
}
\description{
  Compute recursively the ML estimators of location and
  scatter of a multivariate Student t distribution with
  given degrees of freedom, as described in A. Meucci,
  "Risk and Asset Allocation", Springer, 2005, section 4.3
  - "Maximum likelihood estimators"

  Returns \deqn{ \widehat{\mu} = \sum\limits_{t= 1}^{T}
  \frac{w_{t}}{\sum\limits_{s= 1}^{T}w_{s}} x_{t} ,}
  \deqn{\widehat{\Sigma} = \frac{1}{T}\sum\limits_{t=
  1}^{T}w_{t}(x_{t}-\widehat{\mu})(x_{t}-\widehat{\mu})^\prime}

  where, adapted to the Sudent T distribution, \deqn{
  w_{t}\equiv \frac{\nu+N}{\nu+(x_{t}-\widehat{\mu})^\prime
  \widehat{\Sigma}^{-1}( x-\widehat{\mu}) }}
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}, "E 188 -
  Maximum likelihood estimation of a multivariate Student t
  distribution".

  See Meucci's script for "MleRecursionForStudentT.m"
}

