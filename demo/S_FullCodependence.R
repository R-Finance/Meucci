#' This script illustrate the concept of co-dependence, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_FullCodependence.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#############################################################################################################
### Generate draws 
J    = 10000;
N    = 10;
sig2 = 1;

U = runif( J );
X = matrix( NaN, J, N );

for( n in 1 : N )
{
    a = n / 2;
    b = 2 * sig2;
    X[ , n ] = qgamma( U, a, b );
}

NumBins = round( 10 * log( J ));

par( mfrow = c( 3, 1) );
hist( X[ , 1 ], NumBins, xlab = "X_1", main = "histogram of X_1" );
plot( X[ , 1 ], X[ , 2 ], xlab = "X_1", ylab = "X_2" );
hist( X[ , 2 ], NumBins, xlab = "X_2", main = "histogram of X_1" );

