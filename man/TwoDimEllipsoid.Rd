\name{TwoDimEllipsoid}
\alias{TwoDimEllipsoid}
\title{This script computes the location-dispersion ellipsoid of the normalized (unit variance, zero expectation)
first diagonal and off-diagonal elements of a 2x2 Wishart distribution as a function of the inputs,
as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.}
\usage{
  TwoDimEllipsoid(Location, Square_Dispersion, Scale = 1,
    PlotEigVectors = FALSE, PlotSquare = FALSE)
}
\arguments{
  \item{Location}{: [vector] (2 x 1) location vector
  (typically the expected value}

  \item{Square_Dispersion}{: [matrix] (2 x 2) scatter
  matrix Square_Dispersion (typically the covariance
  matrix)}

  \item{Scale}{: [scalar] a scalar Scale, that specifies
  the scale (radius) of the ellipsoid}

  \item{PlotEigVectors}{: [boolean] true then the
  eigenvectors (=principal axes) are plotted}

  \item{PlotSquare}{: [boolean] true then the enshrouding
  box is plotted. If Square_Dispersion is the covariance}
}
\value{
  E : [figure handle]
}
\description{
  This script computes the location-dispersion ellipsoid of
  the normalized (unit variance, zero expectation) first
  diagonal and off-diagonal elements of a 2x2 Wishart
  distribution as a function of the inputs, as described in
  A. Meucci, "Risk and Asset Allocation", Springer, 2005,
  Chapter 2.
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://} See Meucci's script for "TwoDimEllipsoid.m"
}

