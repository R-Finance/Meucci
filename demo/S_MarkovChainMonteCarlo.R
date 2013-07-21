#' This script illustrates the Metropolis-Hastings algorithm, as described in A. Meucci,"Risk and Asset Allocation",
#' Springer, 2005,  Chapter 7.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_MarkovChainMonteCarlo.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Set-up target and candidate
# kernel of the target distribution
kernel = function(x) dnorm( x, 1, 0.5 );

# parameters of the normal candidate distribution
mu  = 0;
sig = 5;

##################################################################################################################
### Set up MH algorithm
nSim = 10000;
xt = matrix( NaN, nSim, 1);
xt[ 1 ] = 0;
nacc = 0;
for( i in 2 : nSim )
{
    # normal candidate
    r = mu + sig * rnorm(1);
    # kernel at candidate
    f1 = kernel( r );
    # kernel at past
    f2 = kernel( xt[ i-1 ] );
    prob = f1 / f2;
    xt[ i ] = xt[ i-1 ];
    if( prob > 1 || runif(1) > (1 - prob) )
    {
        xt[ i ] = r; 
        nacc = nacc + 1;   
    }  
}
##################################################################################################################
### Post-process output
# acceptance rate
print( nacc / nSim );

# display MCMC over time
dev.new();
plot( xt, type = "l" );

# distribution
dev.new();
hist( xt, log(10*nSim) );
