#' This script simulates univariate lognormal variables, as described in  
#' A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 1.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, 
#' E 25- Simulation of a lognormal random variable".
#'
#' See Meucci's script for "S_LognormalSample.m".
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

##################################################################################################################
### Input parameters

nSim = 10000; 
ExpX = 3;
VarX = 5;

##################################################################################################################
### Generate lognormal sample with above parameters

LP = LognormalMoments2Parameters( ExpX, VarX );
sigma = sqrt( LP$sigma_square );

X = rlnorm( nSim, LP$mu, sigma );

##################################################################################################################
### Plots

# plot over time
dev.new();
plot( X, main = "lognormal sample vs observation time" );

# plot histogram
dev.new();
NumBins = round( 10 * log( nSim ) );
hist( X, NumBins, main = "histogram of lognormal sample" );

# plot empirical cdf
dev.new();
f = ecdf( X );
plot( f, col = "red", main = "cdf of lognormal distribution" );

# plot exact cdf
F = plnorm( 1:10000, LP$mu, sigma );
lines ( 1:10000, F, col = "blue" );
legend( "bottomright", 1.9, c("empirical", "exact"), col = c("red", "blue"), lty = 1, bg = "gray90" );


##################################################################################################################
# plot empirical quantile
dev.new();
u= seq( 0.01, 0.99, 0.01 ); # range of quantiles (values between zero and one)
q = quantile( X, u );
plot( u, q, type = "l", xlab="Grade", ylab="Quantile",  lty = 1, col = "red",  main = "quantile of lognormal distribution" );

# plot exact quantile
Q = qlnorm( u, LP$mu, sigma );
lines( u, Q, type = "l", lty = 1, col = "blue" );
legend( "bottomright", 1.9, c( "empirical", "exact" ), col = c( "red", "blue" ), lty = 1, bg = "gray90" );
