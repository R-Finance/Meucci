library(scatterplot3d);

#' This script script shows that the pdf of the r-th order statistics of a tudent t random variable,
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_OrderStatisticsPdfStudentT.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#################################################################################################################
### Input  
mu = 0;
s  = 1;
nu = 10;
T  = 70;

#################################################################################################################
### Pdf of r-th order statistic concentrated around the r/T quantile

rs = 1: T;
x = mu + s * seq( -4, 4, 0.01);

F = pt((x - mu) / s, nu);
f = 1 / s * dt((x - mu) / s, nu);

a = scatterplot3d( 0,  0 , 0, xlim = c(-4 , 4 ), ylim = c( 0, 1 ), zlim = c( 0, 3), xlab = "x", ylab = "r/T", zlab = "pdf" );

for ( n in 1 : length( rs ) )
{
    r = rs[ n ];    
    pdf_rT = gamma( T + 1 ) / ( gamma( r ) * gamma( T - r + 1 )) * ( F ^ (r - 1) ) * (( 1 - F ) ^ ( T - r) ) * f;
    q = mu + s * qt( r / T, nu );
    a$points3d( x, r / T + 0 * x, pdf_rT, type = "l" );
    a$points3d( q, r / T, 0 );
}