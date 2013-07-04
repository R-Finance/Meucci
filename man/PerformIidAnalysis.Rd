\name{PerformIidAnalysis}
\alias{PerformIidAnalysis}
\title{This function performs simple invariance (i.i.d.) tests on a time series, as described in
A. Meucci "Risk and Asset Allocation", Springer, 2005}
\usage{
  PerformIidAnalysis(Dates = dim(Data, 1), Data, Str = "")
}
\arguments{
  \item{Dates}{: [vector] (T x 1) dates}

  \item{Data}{: [matrix] (T x N) data}

  \item{Starting_Prices}{: [vector] (N x 1)}
}
\description{
  This function performs simple invariance (i.i.d.) tests
  on a time series, as described in A. Meucci "Risk and
  Asset Allocation", Springer, 2005
}
\note{
  it checks the evolution over time
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "PerformIidAnalysis.m"
}

