#' This script simulate univariate normal variables, as described in  
#' A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 1.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_NormalSample.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

##################################################################################################################
### Input parameters
nSim   = 10000;
mu     = 3;
sigma2 = 5;

##################################################################################################################
### Generate normal sample with above parameters
sigma = sqrt( sigma2 );
X = rnorm( nSim, mu, sigma);

##################################################################################################################
### Plot the sample# plot over time
dev.new();
plot( X, main = "normal sample vs observation time" );


##################################################################################################################
### Plot the histogram
dev.new();
NumBins = round( 10 * log( nSim ) );
hist( X, NumBins, main = "histogram of normal sample" );

##################################################################################################################
### Compare empirical with exact cdfs

# plot empirical cdf
dev.new();

f = ecdf( X );
plot( f, col = "red", main = "cdf of normal distribution" );

# plot exact cdf
F = pnorm( 1:10000, mu, sigma );
lines ( 1:10000, F, col = "blue" );
legend( "bottomright", 1.9, c("empirical", "exact"), col = c("red", "blue"), lty = 1, bg = "gray90" );

##################################################################################################################
### Compare empirical and exact quantiles

# plot empirical quantile
dev.new();
u= seq( 0.01, 0.99, 0.01 ); # range of quantiles (values between zero and one)
q = quantile( X, u );
plot( u, q, type = "l", xlab="Grade", ylab="Quantile",  lty = 1, col = "red",  main = "quantile of normal distribution" );

# plot exact quantile
Q = qnorm( u, mu, sigma );
lines( u, Q, type = "l", lty = 1, col = "blue" );
legend( "bottomright", 1.9, c( "empirical", "exact" ), col = c( "red", "blue" ), lty = 1, bg = "gray90" );
