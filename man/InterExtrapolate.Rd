\name{InterExtrapolate}
\alias{InterExtrapolate}
\title{Interpolate and extrapolate using n-linear interpolation (tensor product linear).}
\usage{
  InterExtrapolate(V, Xi, nodelist = NULL, method = NULL)
}
\arguments{
  \item{V}{[array] p-dimensional array to be
  interpolated/extrapolated at the list of points in the
  array Xi.}

  \item{Xi}{[array] (n x p) array of n points to
  interpolate/extrapolate. Each point is one row of the
  array Xi.}

  \item{nodelist}{[cell array] (optional) cell array of
  nodes in each dimension.}

  \item{method}{[string] (optional) chacter string, denotes
  the interpolation method used. default method = 'linear'}
}
\value{
  Vpred [array] (n x 1) array of interpolated/extrapolated
  values
}
\description{
  Interpolate and extrapolate using n-linear interpolation
  (tensor product linear).
}
\note{
  Initially written by John D'Errico.

  Extrapolating long distances outside the support of V is
  rarely advisable.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}.

  See Meucci's script for "InterExtrapolate.R"
}

