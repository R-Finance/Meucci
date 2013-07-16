library(MASS);
library(Matrix);
#' This script computes the expected shortfall and the contributions to ES from each factor in simulations, using 
#' the conditional expectation definition of the contributions as described in A. Meucci,"Risk and Asset Allocation",
#' Springer, 2005,  Chapter 5. 
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_ESContributionFactors.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

###################################################################################################################
### Inputs 

N = 30; # number of securities
K = 10; # number of factors
a = runif(N); # allocation
c = 0.95;      # ES confidence

###################################################################################################################
### Generate market simulations
# factor loadings 
B = matrix( runif( N*K ), N, K ); 

# student t parameters
nu  = 10;
eps = 0.1;
g   = 0.1;

A = matrix( runif( K^2 ), K, K ) - 0.5;
C = A %*% t(A);
sigma_f = cov2cor(C);
sigma_u = diag( 1, N );
for( n in 1:(N - 1) )
{
    sigma_u = sigma_u + exp( -g * n ) *( rbind(cbind(matrix(0, N-n, n), diag(array( 1, N - n))), matrix(0,n,N)) +
     rbind(matrix( 0, n, N ),cbind( diag( array( 1, N - n )), matrix( 0, N-n, n)) ));
}

sigma = as.matrix(.bdiag(list( eps * sigma_f, eps^2 * sigma_u))) #block diagonal matrix
corr  = cov2cor( sigma );
diag_sigma = sqrt( diag( sigma ) );
# scenarios
nSim = 10000;
l = matrix( 1, nSim);
X = rmvt( nSim/2,  corr, nu );
X = rbind( X, -X ); # symmetrize simulations
X = X %*% diag( diag_sigma );
X = exp( X );
F = X[ , 1:K ];
U = X[ , (K+1):ncol(X)  ]; 
U = U - l %*% apply( U, 2, mean );
M = F %*% t( B ) + U;

###################################################################################################################
### Risk management
# compute the objective
Psi = M %*% a; 

# compute ES
th = ceiling((1-c) * nSim); # threshold
spc = matrix( 0, nSim, 1 );
spc[ 1 : th ] = 1;
spc = spc / sum( spc );

Sort_Psi = sort( Psi );
Index = order( Psi );
ES_simul = t(Sort_Psi) %*% spc;

# augment factor set to include residual
F_ = cbind( F, U%*%a );
# compute portfolio-level loadings
b_ = cbind( t(a)%*%B, 1 );

# sort factors according to order induced by objective's realizations

Sort_F_ = F_[Index, ];
DES_simul = matrix( NaN, 1, K+1 );
for( k in 1 : (K+1) )
{
    # compute gradient as conditional expectation
    DES_simul[ k ] = t(spc) %*% Sort_F_[ , k ];
}
# compute contributions
ContrES_simul = b_ * DES_simul;

###################################################################################################################
### Plots
dev.new();
bar = barplot( ContrES_simul, main = "contributions to ES" );
axis( 1 , at = bar, labels= 1:ncol(ContrES_simul));
