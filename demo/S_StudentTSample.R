#' This script simulate univariate Student-t variables as described in  
#' A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 1.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_StudentTSample.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

##################################################################################################################
### Input parameters

nSim   = 10000; 
sigma2 = 6;
ExpX   = 2;
VarX   = 7;

##################################################################################################################
### Determine mu, nu and sigma

mu = ExpX;
nu = 2 / (1 - sigma2 / VarX );
sigma = sqrt( sigma2 );

##################################################################################################################
### Generate Student t sample with above parameters using built-in generator

X_a = mu + sigma * rt( nSim, nu );  

##################################################################################################################
### Generate Student t sample with above parameters using stochastic representation

Y = rnorm( nSim, 0, sigma );
Z = rchisq( nSim, nu );
X_b = mu + Y / sqrt( Z / nu );

##################################################################################################################
### Generate Student t sample with above parameters using grade inversion

U = runif( nSim );
X_c = mu + sigma * qt( U, nu );

##################################################################################################################
### Plot histograms
NumBins = round(10 * log(nSim));

par( mfrow = c( 3, 1) );

hist( X_a, NumBins, main = "built-in generator" );
hist( X_b, NumBins, main = "stoch. representation" );
hist( X_c, NumBins, main = "grade inversion" );


#axisLimits = [min(axisLimits(:, 1)), max(axisLimits(:, 2)), min(axisLimits(:, 3)), max(axisLimits(:, 4))];
#subplot(3, 1, 1), axis(axisLimits);
#subplot(3, 1, 2), axis(axisLimits);
#subplot(3, 1, 3), axis(axisLimits);

##################################################################################################################
### Compare empirical quantiles of the three simuations
u= seq( 0.01, 0.99, 0.01 ); # range of quantiles (values between zero and one) = 0.01 : 0.01 : 0.99;  # range of quantiles (values between zero and one)
q_a = quantile( X_a, u );
q_b = quantile( X_b, u );
q_c = quantile( X_c, u );

##################################################################################################################
### Superimpose the the plots of the empirical quantiles

plot( u, q_a, type = "l", xlab="Grade", ylab="Quantile",  lty = 1, col = "red", main = "quantile of Student-t distribution" );
lines( u, q_b, type = "l", lty = 1, col = "blue" );
lines( u, q_c, type = "l", lty = 1, col = "green" );
legend( "bottomright", 1.9, c( "built-in generator", "stoch. representation", "grade inversion" ), col = c( "red" , "blue", "green"),
	lty = 1, bg = "gray90" );