#' This script demonstrates the recursive ML estimation of the location and scatter parameters of a multivariate 
#' Student t distribution, as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_FitSwapToStudentT.m"
#'
#' TO DO: Change colors from TwoDimEllipsoid in each iteration
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
data("usSwapRates" );

##################################################################################################################
### Inputs
ChooseRates = c( 1, 2 ); # 1=2yr; 2=5yr; 3=10yr

Y = cbind( UsSwapRates[ , 1 ], UsSwapRates[ , 3 ] );
X = Y[ -1, ] - Y[ -nrow( Y ), ];

Nus = c( 3, 100 );

#################################################################################################################
### Computations
Tolerance = 10^(-10);
Estimate = array( list() , length(Nus));
for( q in 1 : length(Nus) )
{ 
    Estimate[[q]] = MleRecursionForStudentT( X, Nus[ q ], Tolerance );
}

#################################################################################################################
### Figures
dev.new();
h = plot( X[ , 1 ], X[ , 2 ], bg = "blue", xlab = colnames(UsSwapRates)[1], ylab = colnames(UsSwapRates)[3]);

for( q in 1 : length(Nus) )
{
    M = Estimate[[ q ]]$Mu;
    S = Estimate[[ q ]]$Sigma * Nus[ q ] / ( Nus[ q ] - 2);
    dd = TwoDimEllipsoid( M, S, 2, 0, 0 );
    #set(dd, 'color', 0.7*[rand() rand() rand()]);
}


