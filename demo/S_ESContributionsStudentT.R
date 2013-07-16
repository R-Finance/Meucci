library(MASS);
library(Matrix);
#' This script computes the expected shortfall and the contributions to ES from each security: 
#' - analytically, under the Student t assumption for the market
#' - in simulations, using the conditional expectation definition of the contributions
#' Described in A. Meucci,"Risk and Asset Allocation",Springer, 2005,  Chapter 5.  
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_ESContributionsStudentT.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs 
# number of assets
N = 40; 

# market parameters (student t distribution)
Mu = runif(N);
A  = matrix( runif( N^2 ), N, N )-.5;
Sigma = A %*% t(A);
nu = 7;

# allocation
a = runif(N) - 0.5;

# ES confidence
c = 0.95;

nSim = 10000;

###################################################################################################################
### Generate market scenarios
l = matrix( 1, nSim, 1);
diag_s    = diag( sqrt( diag( Sigma ) ) );
invdiag_s = diag( 1 / sqrt( diag( Sigma ) ) );
C = invdiag_s * Sigma * invdiag_s;
X = rmvt( nSim / 2, C, nu);
X = rbind( X, -X ); # symmetrize simulations
M = l %*% t( Mu ) + X %*% diag_s;

###################################################################################################################
### Simulations
# compute the objective
Psi = M %*% a; 

# compute cut-off spectrum (step function) for empirical ES estimation, see (5.218)
th = ceiling((1-c) * nSim); # threshold
spc = matrix( 0, nSim, 1 );
spc[ 1 : th ] = 1;
spc = spc / sum( spc );

# compute ES
Sort_Psi = sort( Psi );
Index = order( Psi );
ES_simul = t(Sort_Psi) %*% spc;

# sort market according to order induced by objective's realizations
Sort_M = M[ Index, ];
DES_simul = matrix(NaN, 1, N);
for( n in 1 : N )
{
    # compute gradient as conditional expectation
    DES_simul[ n ] = t(spc) %*% Sort_M[ , n ];
}

# compute contributions
ContrES_simul = a * t( DES_simul);

###################################################################################################################
### Analytical
# this does NOT depend on the allocation...
ES_standardized = 1 / ( 1 - c ) * integrate(qt, lower=10^(-8) ,upper=(1-c), df=7)$value;

# ...the dependence on the allocation is analytical
ES_an  = t( Mu ) %*% a + ES_standardized * sqrt(t(a) %*% Sigma %*% a);
DES_an = Mu + ES_standardized * Sigma %*% a / sqrt( t(a) %*% Sigma %*% a)[1];
ContrES_an = a * DES_an;

###################################################################################################################
### Plots
dev.new();
par( mfrow = c( 2, 1 ) );
bar = barplot(t(ContrES_an),  main = "contributions to ES, analytical" );
axis( 1 , at = bar, labels= 1:length(ContrES_an[,1]));

bar = barplot(t(ContrES_simul), main = "contributions to ES, simulations" );
axis( 1 , at = bar, labels= 1:length(ContrES_simul[,1]));