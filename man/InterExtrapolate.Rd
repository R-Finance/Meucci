\name{InterExtrapolate}
\alias{InterExtrapolate}
\title{Interpolate and extrapolate using n-linear interpolation (tensor product linear).}
\usage{
  InterExtrapolate(V, Xi, nodelist, method, ...)
}
\arguments{
  \item{V}{: [array] p-dimensional array to be
  interpolated/extrapolated at the list of points in the
  array Xi.}

  \item{Xi}{: [array] (n x p) array of n points to
  interpolate/extrapolate. Each point is one row of the
  array Xi.}

  \item{nodelist}{: [cell array] (optional) cell array of
  nodes in each dimension.}

  \item{method}{: [string] (optional) chacter string,
  denotes the interpolation method used. default method =
  'linear'}
}
\value{
  Vpred : [array] (n x 1) array of
  interpolated/extrapolated values
}
\description{
  Interpolate and extrapolate using n-linear interpolation
  (tensor product linear).
}
\note{
  Initially written by John D'Errico Vpred = interpne(V,Xi)
  Vpred = interpne(V,Xi,nodelist) Vpred =
  interpne(V,Xi,nodelist,method) Extrapolating long
  distances outside the support of V is rarely advisable.
}
\examples{
[x1,x2] = meshgrid(0:.2:1);
 z = exp(x1+x2);
 Xi = rand(100,2)*2-.5;
 Zi = interpne(z,Xi,{0:.2:1, 0:.2:1},'linear');
 surf(0:.2:1,0:.2:1,z)
 hold on
 plot3(Xi(:,1),Xi(:,2),Zi,'ro')
}
\author{
  Xavier Valls \email{flamejat@gmail.com}
}
\references{
  \url{http://symmys.com/node/170} See Meucci's script for
  "InterExtrapolate.R"
}

