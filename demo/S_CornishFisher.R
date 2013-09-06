#'This script compares the Cornish-Fisher estimate of the VaR with the true analytical VaR under the lognormal 
#'assumptions as described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 5.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_CornishFisher.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

###################################################################################################################
### Inputs 
mu  = 0.05; 
sig = 0.05; # NB: change here and see the impact of approximation

###################################################################################################################
### Process data
E_X  = exp( mu + sig ^ 2 / 2 );
Sd_X = exp( mu + sig ^ 2 / 2 ) * sqrt( exp( sig ^ 2 ) - 1 );
Sk_X = sqrt( exp( sig ^ 2 ) - 1 ) * ( exp( sig ^ 2 ) + 2 );

c = seq(0.001, 0.999, 0.001 );
z = qnorm( c );

Q_CF = E_X + Sd_X * ( z + Sk_X  /  6 * ( z ^ 2 - 1 ) );
Q_true = qlnorm( c,mu,sig );

x = Q_true;
f = dlnorm( x, mu, sig );

###################################################################################################################
### Plots
dev.new();
plot( x, f, type= "l", main = "pdf" );

dev.new();
plot( c, Q_true, type = "l", col = "red", main = "quantile" );
lines( c, Q_CF );
legend( "topleft", 1.9, c( "true", "Cornish-Fisher approx" ), col = c( "black","red" ),
     lty=1, bg = "gray90" );
