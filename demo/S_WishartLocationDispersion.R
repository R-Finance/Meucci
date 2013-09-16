#' This script computes the location-dispersion ellipsoid of the normalized (unit variance, zero expectation)
#' first diagonal and off-diagonal elements of a 2x2 Wishart distribution as a function of the inputs,
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, 
#' "E 87 - Correlation and location-dispersion ellipsoid", "E 75 - Simulation of a Wishart random variable".
#'
#' See Meucci's script for "S_WishartLocationDispersion.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' 

###################################################################################################################
### Set input parameters
s  = c( 1, 1 ); # variances
r  = -0.9; # correlation
Sigma = diag( s ) %*% rbind( c( 1, r ), c( r, 1 ) ) %*% diag( s );
nu = 5; # degrees of freedom
nSim = 10000;

###################################################################################################################
### Set input parameters

W_xx   = matrix( NaN, nSim, 1 ); 
W_yy   = matrix( NaN, nSim, 1 ); 
W_xy   = matrix( NaN, nSim, 1 ); 
Vec_W  = matrix( NaN, nSim, 4 ); 
Dets   = matrix( NaN, nSim, 1 ); 
Traces = matrix( NaN, nSim, 1 ); 


for( j in 1 : nSim )
{
  X = rmvnorm( nu, matrix( 0, 2, 1 ), Sigma);
  W = t( X ) %*% X;
 
  Dets[ j ] = det( W );
  Traces[ j ] = tr( W );

  W_xx[ j ] = W[ 1, 1 ];
  W_yy[ j ] = W[ 2, 2 ];
  W_xy[ j ] = W[ 1, 2 ];

  Vec_W [ j, ] = as.vector( W );
}

# compute expected values of W_xx and W_xy, see (2.227) in "Risk and Asset Allocation - Springer
E_xx = nu * Sigma[ 1, 1 ];
E_xy = nu * Sigma[ 1, 2 ];

# compute covariance matrix of W_xx and W_xy, see (2.228) in "Risk and Asset Allocation - Springer
m = 1; 
n = 1; 
p = 1; 
q = 1;
var_Wxx = nu * ( Sigma[ m, p ] * Sigma[ n, q ] + Sigma[ m, q ] * Sigma[ n, p ] );
m = 1;
n = 2;
p = 1;
q = 2;
var_Wxy = nu * ( Sigma[ m, p ] * Sigma[ n, q ] + Sigma[ m, q ] * Sigma[ n, p ] );
m = 1; 
n = 1;
p = 1;
q = 2;
cov_Wxx_Wxy = nu * ( Sigma[ m, p ] * Sigma[ n, q ] + Sigma[ m, q ] * Sigma[ n, p ] );

S_xx_xy = rbind( cbind( var_Wxx, cov_Wxx_Wxy ), cbind( cov_Wxx_Wxy, var_Wxy ));

# compute X_1 and X_2, i.e. normalized version of W_xx and W_xy
X_1 = ( W_xx - E_xx ) / sqrt( var_Wxx );
X_2 = ( W_xy - E_xy ) / sqrt( var_Wxy );
X = cbind( X_1, X_2 );

# compute expected value and covariance of X_1 and X_2
E = rbind( 0, 0 );
E_hat = t( apply( X, 2, mean) );

S = diag( 1 / c( sqrt( var_Wxx ), sqrt( var_Wxy ))) %*% S_xx_xy %*% diag( 1 / c( sqrt( var_Wxx ), sqrt( var_Wxy )));
S_hat = cov( X );

dev.new();
plot( X_1, X_2, xlab = "X_1", ylab = "X_2");

TwoDimEllipsoid(E, S, 1, TRUE, FALSE);
