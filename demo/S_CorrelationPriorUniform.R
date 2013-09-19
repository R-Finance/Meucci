#' This script shows how a jointly uniform prior on the correlations implies that the marginal distribution of 
#' each correlation is peaked around zero , as described in A. Meucci,"Risk and Asset Allocation",Springer, 2005, 
#'  Chapter 7.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 281 - Bayesian: prior on correlation".
#'
#' See Meucci's script for "S_CorrelationPriorUniform.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs
N = 3; # dimensionality of the problem
K = N * (N - 1) / 2;
J = 10000; # number of simulations

##################################################################################################################
### Compute correlations in all scenarios
CorrsAsTensor = array(0, dim = c(J,N,N) );
Eigs = NULL;
j    = 1;

while( j < J )
{
    C    = 2 * matrix( runif(K), 1, K ) - 1;
    Corr = diag( 1, N );
    k    = 0;
    for( n in 1 : ( N - 1 ) )
    {
        for( m in ( n + 1 ) : N )
        {
            k = k + 1;
            Corr[ n, m ] = C[ k ];
            Corr[ m, n ] = Corr[ n, m ];
        }
    }
    E = eigen(Corr)$values;
    if( min(E) > 0 )
    {
        CorrsAsTensor[ j, , ] = Corr;
        j = j + 1;
    }
    Eigs = rbind(Eigs, t(E) );
}

#####################################################################################################################
### Reassemble results in an entry-wise structure that runs on the upper triangular portion of the correlation
CorrsAsEntries = NULL;
for( n in 1 : (N-1) )
{
    for( m in (n + 1) : N )
    {
        el = list( Values = CorrsAsTensor[ , n, m ], Names =  paste("{", n, m,"}"));
              #  CorrsAsEntries[ k ]$Names  = ["\theta_{") num2str(n) num2str(m) "}"];
        CorrsAsEntries = rbind( CorrsAsEntries, el )
    }
}
#####################################################################################################################
### Plots
# univariate marginals
K     = nrow( CorrsAsEntries );
Nbins = round( 5 * log( J ) );
for( k in 1 : K )
{
    dev.new();
    hist(CorrsAsEntries[ k, ]$Values, Nbins, xlab = "", ylab = "", main = bquote( paste("histogram of ", theta, .(CorrsAsEntries[k,]$Names))));
}

# bivariate marginals
for( k in 1 : (K-1) )
{
    for( j in (k + 1) : K )
    {
        dev.new();
        plot(CorrsAsEntries[k, ]$Values,CorrsAsEntries[j, ]$Values, xlab = "", ylab = "",
         main = paste( CorrsAsEntries[ k ]$Names, ' - ', CorrsAsEntries[ j ]$Names ));
    }
}
