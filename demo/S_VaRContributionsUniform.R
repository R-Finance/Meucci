
#' This script computes the VaR and the contributions to VaR from each security 
#'  - analytically, under the elliptical-uniform assumption for the market
#'  - in simulations, using the conditional expectation definition of the contributions
#' Described in A. Meucci,"Risk and Asset 
#' Allocation",Springer, 2005,  Chapter 5.  
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_VaRContributionsUniform.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

###################################################################################################################
### Inputs 
# number of assets
N = 10; 

# market parameters (uniform on ellipsoid)
Mu = matrix(runif(N));
A  = matrix( runif(N*N), N, N) - 0.5;
Sigma = A * t(A);

# allocation
a = matrix( runif(N) ) - 0.5;

# quantile level
c = 0.95;

nSim = 10000;

###################################################################################################################
### Generate market scenarios
X = GenerateUniformDrawsOnUnitSphere(nSim, N); # uniform on sphere
M = matrix( 1, nSim, 1) %*% t(Mu) + X %*% t(A);

###################################################################################################################
### Compute contributions by simulations (brute-force approach)
# compute and sort the objective
Psi = M %*% a;
Q_sim = quantile( Psi, (1 - c) );

e = mean( abs( a )) / 100; # perturbation
DQ_simul = matrix( NaN, 1, N) ;
for( n in 1 : N )
{
    # compute gradient
    a_e = a;
    a_e[ n ] = a[ n ] + e;
    
    Psi_e = M %*% a_e;
    Q_sim_e = quantile(Psi_e, (1 - c) );
    DQ_simul[ n ] = ( Q_sim_e - Q_sim )/e;
}
# compute contributions
ContrQ_simul = a * t( DQ_simul );

###################################################################################################################
### Compute contributions analytically
# compute quantile of standardized marginal (1-dim generator) in simulations this does NOT depend on the allocation...
gc = quantile(X[  ,1 ], (1 - c));

# ...the dependence on the allocation is analytical
Q_an  = t(Mu) %*% a + gc * sqrt( t(a) %*% Sigma %*% a );
DQ_an = Mu + gc * Sigma %*% a / sqrt( t(a) %*% Sigma %*% a )[1];
ContrQ_an = a * DQ_an;

###################################################################################################################
# plots
dev.new();
par( mfrow = c(2,1) );
bar = barplot( t(ContrQ_an), xlim = c(0, N+1), main = "contributions to VaR, analytical" );
axis( 1 , at = bar, labels= 1:(N+1));

bar = barplot( t(ContrQ_simul), xlim = c(0, N+1), main = "contributions to VaR, simulations" );
axis( 1 , at = bar, labels= 1:(N+1));
