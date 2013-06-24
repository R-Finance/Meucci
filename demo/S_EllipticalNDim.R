library(mvtnorm);
#'This script decomposes the N-variate normal distribution into its radial and uniform components
#' then it uses the uniform component to generate an elliptical distribution with location parameter 
#' Mu and dispersion parameter Sigma, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_EllipticalNDim.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

##################################################################################################################
### Parameters

N = 30;
nSim = 10000;
nu = 0.1;
t2 = 0.04;

##################################################################################################################
### Random matrix and elliptical draws
Mu = runif( N );
A  = matrix( runif( N * N ), c( N, N )) - 0.5;
Sigma = A %*% t( A );

Y = rmvnorm( nSim, matrix( 0, N, 1 ), diag( 1, N ));

# radial distribution (normal case ~ square root of chi-square with N degrees of freedom)
R = matrix( sqrt( apply( Y * Y, 1, sum )));


# uniform distribution on unit sphere
U = Y / ( R %*% matrix( 1, 1, N )); 

tau = sqrt( t2 );
R_New = rlnorm( nSim, nu, tau );

# N-variate elliptical distribution 
X = matrix( 1, nSim, 1 ) %*% t( Mu ) + ( R_New %*% matrix( 1, 1, N )) * ( U %*% t( A )); 

##################################################################################################################
### Plots 
# visualize projection on m-n coordinates

m = 1;
n = 3;
xlabel = paste( "X_" , m );
ylabel = paste( "X_", n );
dev.new();
plot( X[ , m ], X[ , n ], xlab = xlabel, ylab = ylabel);

# visualize n-th marginal
n = 4;
xlabel = paste( "X_", m );
NumBins = round(10 * log(nSim));
dev.new();
hist( X[ , n ], NumBins, xlab = xlabel, main= "histogram");

