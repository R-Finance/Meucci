#' This file generates paths for a volatility clustering, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_VolatilityClustering.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
##################################################################################################################
### Input parameters
mu = 0.05; # mean
a  = 0.03;
b  = 0.96;
s  = 0.01;
T  = 1000;

##################################################################################################################
### Simulate path
z = rnorm(T);
s2 = s^2;
eps = array( NaN, T );
eps[ 1 ] = s2;
for( t in 1 : (T - 1) )
{
    s2[ t + 1 ]  = s^2 + a * ( z[ t ]^2) + b * s2[ t ];
    eps[ t + 1 ] = mu + sqrt( s2[ t + 1] ) * z[ t + 1 ];
}

dev.new();
plot(eps, type = "l", main = "GARCH(1,1) process vs. time", xlab = "", ylab = "" );
