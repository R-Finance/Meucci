#' This script computes the quantile (VaR) : 
#' - analytically, under the Student t assumption for the market.
#' - in simulations, using the sample quantile.
#' - using the extreme value theory approximation
#' Described in A. Meucci,"Risk and Asset Allocation",Springer, 2005,  Chapter 5.  
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_ExtremeValueTheory.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

if ( !require( "fExtremes" ) ) stop("fExtremes package installation required for this script")
##################################################################################################################
### Market parameters (student t distribution)
m  = 1;
s  = 2;
nu = 7;
nSim = 10000;

th = 0.95; # EVT threshold
c = seq( th , 0.999, 0.001 ); # confidence range for quantiles

###################################################################################################################
### Analytical
Q_an = m + s * qt(1 - c, nu);

###################################################################################################################
### Simulations
# generate objective's scenarios
X = rt( nSim/2, nu);
X = rbind( X, -X );  # symmetrize simulations
Psi = m + s * X;
Q_simul = quantile( Psi, (1 - c));

###################################################################################################################
### EVT approximation
psi_hat = quantile(Psi, (1 - th));
Excess = psi_hat - Psi[ Psi < psi_hat ];
xi_v = gpdFit(Excess);
xi   = xi_v@fit$par.ests[1];
v    = xi_v@fit$par.ests[2];


Fpsi_hat = 1 - th;
Q_EVT = psi_hat + v / xi * ( 1 - ( ( 1 - c ) / Fpsi_hat ) ^ ( -xi ) );

###################################################################################################################
### Plots
dev.new();
plot(c, Q_an, type = "l", xlab = "confidence, c", ylab = "quantile based satisfaction, Q_c(\alpha)" );
lines(c, Q_simul, col = "green" );
lines(c, Q_EVT, col = "red" );
legend( "bottomleft", 1.9, c( "exact", "simulations", "EVT" ), col = c( "black","green", "red" ),
     lty=1, bg = "gray90" );
