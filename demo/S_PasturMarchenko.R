#' This script illustrate the Marchenko-Pastur limit of runifom matrix theory, as described in A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_PasturMarchenko.m"
#'

##################################################################################################################
### Inputs
T = 1500;
N = 900;

##################################################################################################################
### Empirical eigenvalues

#X = matrix( runif(T*N), T, N ) ; # normal
#X = (matrix( runif(T*N), T, N ) - 0.5) * sqrt(12); # uniform
X = log(matrix( runif( T*N ), T, N )) + 1; # exponential

Y = ( t(X) %*% X ) / T; # symmetrize and rescale
E = t(eigen(Y)$values);

NumBins = ceiling( 10 * log( length( E )));
h = hist(E, NumBins, 1); 
t_= h$mids;
b = h$counts;
D = t_[ 2 ] - t_[ 1 ];
h = b / (D * N);

##################################################################################################################
### Theoretical eigenvalues
q = N / T;
t_min = ( 1 - sqrt( q )) ^ 2;
t_max = ( 1 + sqrt( q )) ^ 2;
t = seq(t_[ 1 ], t_[length(t_)], (t_[ length(t_) ]- t_[ 1 ])/100 );
a = pmax( t_max - t, 0);
b = pmax( t - t_min, 0);
y = 1 / ( q * 2 * pi * t) * sqrt(a * b);

##################################################################################################################
### Plots
#barplot(t_,h);
plot(t_,h, type="h", lwd=5);
lines(t , y, col = 'red', lwd = 3);
