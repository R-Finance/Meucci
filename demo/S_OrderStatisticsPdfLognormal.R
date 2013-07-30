library(scatterplot3d);

#' This script script shows that the pdf of the r-th order statistics of a lognormal random variable,
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_OrderStatisticsPdfLognormal.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

if ( !require( "scatterplot3d" ) ) stop("scatterplot3d package installation required for this script")

#################################################################################################################
### Input  

mu = 0.2;
s  = 0.25;
T  = 70;

#################################################################################################################
### Pdf of r-th order statistic concentrated around the r/T quantile

rs = 1 : T;
x  = seq( 0 , 2.5 * exp(mu + s * s / 2), 0.01 );

F = plnorm( x, mu, s );
f = dlnorm( x, mu, s );

#matrix to plot

a = scatterplot3d( 0,  0 , 0, xlim=c(0,4), ylim=c(0,1), zlim=c(0,10), xlab = "x", ylab = "r/T", zlab = "pdf" );

for ( n in 1 : length( rs ) )
{
    r = rs[ n ];    
    pdf_rT = gamma( T + 1 ) / ( gamma( r ) * gamma( T - r + 1 )) * ( F ^ (r - 1) ) * (( 1 - F ) ^ ( T - r) ) * f;
    q = qlnorm( r / T, mu, s );
    a$points3d( x, r / T + 0 * x, pdf_rT );
    a$points3d( q, r / T, 0 );
}