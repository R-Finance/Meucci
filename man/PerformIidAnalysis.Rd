\name{PerformIidAnalysis}
\alias{PerformIidAnalysis}
\title{Performs simple invariance (i.i.d.) tests on a time series.}
\usage{
  PerformIidAnalysis(Dates = dim(Data)[1], Data, Str = "")
}
\arguments{
  \item{Dates}{: [vector] (T x 1) dates}

  \item{Data}{: [matrix] (T x N) data}

  \item{Str}{: [string] title for the plot}
}
\description{
  This function performs simple invariance (i.i.d.) tests
  on a time series, as described in A. Meucci "Risk and
  Asset Allocation", Springer, 2005
}
\note{
  it checks the evolution over time

  it checks that the variables are identically distributed
  by looking at the histogram of two subsamples

  it checks that the variables are independent by looking
  at the 1-lag scatter plot

  under i.i.d. the location-dispersion ellipsoid should be
  a circle
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  A. Meucci - "Exercises in Advanced Risk and Portfolio
  Management" \url{http://symmys.com/node/170}.

  See Meucci's script for "PerformIidAnalysis.m"
}

