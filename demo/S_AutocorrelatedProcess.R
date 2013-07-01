
#' This script simulates a Ornstein-Uhlenbeck AR(1) process, as described in A. Meucci, "
#' Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_AutocorrelatedProcess.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Input parameters
theta = 0.1;  # reversion speed
m     = 0.05; # long term mean 
sigma = 0.01; # volatility
T     = 10^4; # number of steps
tau   = 0.01; # discrete time interval

##################################################################################################################
### Determine parameters
var = sigma^2 / 2 / theta * ( 1 - exp( -2 * theta * tau ) );
sd = sqrt(var);
eps = rnorm( T, 0, sd );

x = matrix( NaN, T, 1);
x[ 1 ] = 0;

for( t in 1 : (T - 1) )
{
    x[ t + 1 ] = m + exp( -theta * tau ) * ( x[ t ] - m ) + eps[ t ];
}

dev.new()
plot( x, type="l", main = "AR(1) process vs. time" );
