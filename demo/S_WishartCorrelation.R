#' This script computes the correlation of the first diagonal and off-diagonal elements 
#' of a 2x2 Wishart distribution as a function of the inputs, as described in A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, 
#' "E 87 - Correlation and location-dispersion ellipsoid", "E 75 - Simulation of a Wishart random variable".
#'
#' See Meucci's script for "S_WishartCorrelation.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' 

###################################################################################################################
### Inputs

s = c( 1, 1);
nu = 15;
rhos = seq( -0.99, 0.99, 0.01 );
nrhos = length(rhos);

###################################################################################################################
### Compute the correlation using simulation

corrs2 = matrix( NaN, nrhos, 1);
for( uu in 1 : nrhos )
{
    rho = rhos[ uu ];
    Sigma = diag( s ) %*% rbind( c( 1, rho ), c( rho, 1 ) ) %*% diag( s );

    # compute expected values of W_xx and W_xy, see (2.227) in "Risk and Asset Allocation - Springer
    E_xx = nu * Sigma[ 1, 1 ];
    E_xy = nu * Sigma[ 1, 2 ];

    # compute covariance matrix of W_xx and W_xy, see (2.228) in "Risk and Asset Allocation - Springer
    m = 1; n = 1; p = 1; q = 1;
    var_Wxx = nu * ( Sigma[ m, p ] * Sigma[ n, q ] + Sigma[ m, q ] * Sigma[ n, p ] );
    
    m = 1; n = 2; p = 1; q = 2;
    var_Wxy = nu * ( Sigma[ m, p ] * Sigma[ n, q ] + Sigma[ m, q ] * Sigma[ n, p ] );
    
    m = 1; n = 1; p = 1; q = 2;
    cov_Wxx_Wxy = nu * ( Sigma[ m, p ] * Sigma[ n, q ] + Sigma[ m, q ] * Sigma[ n, p ] );

    S_xx_xy = rbind( cbind( var_Wxx, cov_Wxx_Wxy ), cbind( cov_Wxx_Wxy, var_Wxy ));

    # compute covariance of X_1 and X_2
    S = diag( 1 / c( sqrt( var_Wxx ), sqrt( var_Wxy ))) %*% S_xx_xy %*% diag( 1 / c( sqrt( var_Wxx ), sqrt( var_Wxy )));

    # correlation = covariance
    corrs2[ uu ] = S[ 1, 2 ];
}

###################################################################################################################
### Analytical correlation
corrs = sqrt( 2 ) * rhos / sqrt( 1 + rhos ^ 2);

dev.new();
plot(rhos, corrs, xlab = expression( paste("input ", rho)), ylab = "Wishart correlation");
lines( rhos, corrs2, col = "red" );
