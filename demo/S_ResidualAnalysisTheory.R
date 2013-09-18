#' This script performs the analysis of residuals, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 117 - Time series factors: analysis of residuals II".
#'
#' See Meucci's script for "S_ResidualAnalysisTheory.m"
#' @note See #' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 116 - Time series factors: analysis of residuals I".
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' 

##################################################################################################################
### Inputs
N    = 3; # number of stocks
K    = 2; # number of factors
tau  = 12 / 12; # years
nSim = 10000;
useNonRandomFactor = FALSE; # NB: change here to true/false
useZeroMeanLinearFactors = TRUE; # NB: change here to true/false

##################################################################################################################
### Generate covariance matrix with volatility 25#
Diag = 0.25 * array( 1, N+K );
# make one factor deterministic -> add constant
if( useNonRandomFactor )
{
    Diag[ length(Diag) ] = 10^( -26 );
}
dd = matrix( runif( (N + K)^2), N + K, N + K ) - 0.5;
dd = dd %*% t( dd );
C = cov2cor( dd );

covJointXF = diag( Diag, length(Diag) ) %*% C %*% diag( Diag, length(Diag) );

SigmaX  = covJointXF[ 1:N, 1:N ];
SigmaXF = covJointXF[ 1:N, (N+1):(N+K) ];
SigmaF  = covJointXF[ (N+1):(N+K), (N+1):(N+K) ];

##################################################################################################################
### Generate mean returns for stock and factors ~ 10# per annum
muX = 0.1 * rnorm(N); 
muF = 0.1 * rnorm(K); 

##################################################################################################################
### Statitics of linear returns analytically, since Y = 1+[R; Z] is lognormal
mu   = c( muX, muF );
E_Y  = exp( mu * tau + diag( covJointXF * tau ) / 2);
E_YY = ( E_Y %*% t(E_Y) ) * exp( covJointXF * tau );
E_R  = E_Y[ 1:N ] - array( 1, N );
E_Z  = E_Y[ (N+1):length(E_Y) ] - array( 1, K );
E_RR = E_YY[ 1:N, 1:N ] - array( 1, N ) %*% t(E_R) - E_R %*% t( array( 1, N ) ) - matrix( 1, N, N );
E_ZZ = E_YY[ (N+1):nrow(E_YY), (N+1):ncol(E_YY) ] - array( 1, K ) %*% t(E_Z) - E_Z %*% t(array( 1, K )) - matrix( 1, K, K );
E_RZ = E_YY[ 1:N, (N+1):ncol(E_YY) ] - array( 1, N ) %*% t(E_Z) - E_R %*% t(array(1,K)) - matrix( 1, N, K );
SigmaZ  = E_ZZ - E_Z %*% t(E_Z);
SigmaR  = E_RR - E_R %*% t(E_R);
SigmaRZ = E_RZ - E_R %*% t(E_Z);

##################################################################################################################
### Generate Monte Carlo simulations
sims = rmvnorm( nSim, mu * tau, covJointXF * tau );
X = sims[ , 1:N ];
F = sims[ , (N+1):ncol(sims) ];
R = exp(X) - 1;
Z = exp(F) - 1;
# enforce Z sample to be zero-mean, equivalent to muF = -diag(SigmaF)/2
if( useZeroMeanLinearFactors )
{
    Z = Z - repmat( apply( Z, 2, mean ), nSim, 1 );
}

##################################################################################################################
### Compute sample estimates
E_R_hat = matrix( apply( R, 2, mean) );
E_Z_hat = matrix( apply( Z, 2, mean) );
SigmaR_hat = ( dim(R)[1] - 1 ) / dim(R)[1] * cov( R );
SigmaZ_hat = ( dim(Z)[1] - 1 ) / dim(Z)[1] * cov( Z );

##################################################################################################################
### Compute simulation errors
errSigmaR = norm( SigmaR - SigmaR_hat, "F" ) / norm( SigmaR, "F" );
printf( "Simulation error in sample cov(R) as a percentage on true cov(R) = %0.1f \n", errSigmaR * 100 );
errSigmaZ = norm( SigmaZ - SigmaZ_hat, "F" ) / norm( SigmaZ, "F" );
printf( "Simulation error in sample cov(Z) as a percentage on true cov(Z) = %0.1f \n", errSigmaZ * 100 );

##################################################################################################################
### Compute OLS loadings for the linear return model
B = E_RZ %*% solve( E_ZZ );
ZZ = t(Z) %*% Z;
B_hat = t(R) %*% Z %*% solve(ZZ);
errB = norm( B - B_hat, "F" ) / norm( B, "F" );
printf( "Simulation error in sample OLS loadings as a percentage on true OLS loadings = %0.1f \n", errB * 100 );

U = R - Z %*% t( B_hat );
Corr = cor(cbind( U, Z ));

Corr_U  = Corr[ 1:N, 1:N ];
Corr_UZ = Corr[ 1:N, (N+1):(N+K) ];

SigmaU_hat = ( dim(U)[1] - 1 ) / dim(U)[1] * cov( U );
BSBplusSu  = B_hat %*% SigmaZ_hat %*% t( B_hat ) + SigmaU_hat;
errSigmaR_model1 = norm( SigmaR_hat - BSBplusSu, "F" ) / norm( SigmaR_hat, "F" );
printf( "Sigma_R -( B * Sigma_Z * t(B) + Sigma_U) = %0.1f \n", errSigmaR_model1 * 100 );

