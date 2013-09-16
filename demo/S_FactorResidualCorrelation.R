#' This script illustrates exogenous loadings and endogenous factors the true analytical VaR under the lognormal  
#' assumptions from the estimation interval to the investment horizon, as described in A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_FactorResidualCorrelation.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Input parameters:
N     = 4; # market size
nSim  = 10000;
mu    = 0.1 + 0.3 * runif(N);
sigma = 0.5 * mu; 
dd = matrix(rnorm( N*N ), N, N );
Corr = cov2cor( dd %*% t( dd ) );
Sigma = diag( sigma, length(sigma) ) %*% Corr %*% diag( sigma, length(sigma) ); 

##################################################################################################################
### Generate simulations for X
X = MvnRnd(mu, Sigma, nSim);

##################################################################################################################
### Generate a random vector beta
beta = matrix(1, N ) + rnorm(N) * 0.1;

##################################################################################################################
### Compute factor realization by cross-sectional regression and residuals
F = ( X %*% beta ) / ( t( beta ) %*% (beta) )[1];
F = solve( t( beta ) %*% beta)[1] * ( X %*% beta );

# compute residual
U = X - F %*% t(beta);

# correlation of residuals U among themselves and with factors F
R = cor( cbind( F, U ) );
print(R);

