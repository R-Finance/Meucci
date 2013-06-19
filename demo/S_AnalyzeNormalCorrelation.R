#' This script considers a bivariate normal market and display the correlation and the condition number of the
#' covariance matrix, as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_AnalyzeNormalCorrelation.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

###################################################################################################################
### Set input parameters

Mu = rbind( 0, 0 )
s  = c( 1, 1 );

rhos = seq( -0.99, 0.99, 0.01 );
nrhos = length( rhos );

Cs = array( NaN, nrhos ); 
CRs = array( NaN, nrhos ); 


###################################################################################################################
### Iterate of rho values and compute the correlation and condition numberfor ( n in 1 : nrhos )

for ( n in 1 : nrhos )
{
	rho = rhos[ n ] ;
    Sigma = rbind( c(s[1]^2, rho * s[1] * s[2]), c(rho * s[1] * s[2], s[2]^2) );

    Covariance = Sigma;
    Standard_Deviation = sqrt( diag( Covariance ) );
    Correlation = diag( 1 / Standard_Deviation ) %*% Covariance %*% diag( 1 / Standard_Deviation );
    
    Lambda = eigen( Covariance );

    Cs[n]  = Correlation[ 1, 2 ];
    CRs[n] = min( Lambda$values ) / max( Lambda$values );
}

###################################################################################################################
### Display the results
par( mfrow = c( 2, 1) );
plot( rhos, Cs, xlab = "r", ylab = "correlation", type ="l" );
plot( rhos, CRs, xlab = "r", ylab = "condition ratio", type ="l" );
