#' This script generates draws for the sum of random variables, as described in  
#' A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 1.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 22- Sum of random variables via simulation".
#'
#' See Meucci's script for "S_NonAnalytical.m".
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Input parameters

nSim  = 10000;
mu_St = 0;
s2_St = 0.1;
nu_St = 8; # NOTE: see how the final results change if you increase nu (=4 is enough)

mu_LN = 0.1;
s2_LN = 0.2;

##################################################################################################################
### Generate draws

# generate Student sample with above parameters
s_St = sqrt( s2_St );
X = mu_St + s_St * rt( nSim, nu_St );

# generate lognormal sample with above parameters
s_LN = sqrt( s2_LN );
Y = rlnorm( nSim, mu_LN, s_LN );

# sum samples
Z = X + Y;

##################################################################################################################
### Plot the sample Z
dev.new();
plot( Z,  xlab="simulations", ylab="Z", main = "sample vs observation time" );

##################################################################################################################
### Plot the histogram of Z
dev.new();
NumBins = round( 10 * log( nSim ) );
hist( Z, NumBins, xlab="Z", main="sample histogram" );

##################################################################################################################
### Plot the empirical cdf of Z
dev.new();
f = ecdf( Z );
plot( f, xlab="Z", main="empirical cdf" );

##################################################################################################################
### Plot the empirical quantile of Z
dev.new();
u= seq( 0.01, 0.99, 0.01 ); # range of quantiles (values between zero and one)
q = quantile( Z, u );
plot( u, q, type = "l", xlab="Grade", ylab="Quantile",  lty = 1,  main = "empirical quantile" );
