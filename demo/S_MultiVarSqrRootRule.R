#' This script illustrates the multivariate square root rule-of-thumb 
#' Described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 95 - Multivariate square-root rule".
#'
#' See Meucci's script for "S_MultiVarSqrRootRule.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
data("swaps");

##################################################################################################################
### Aggregation steps in days
Steps = c( 1, 5, 22 ); 

##################################################################################################################
### Plots
Agg = list();

dev.new();
plot( swaps$X[ , 1 ],swaps$X[ , 2], xlab = swaps$Names[[1]][1], ylab = swaps$Names[[2]][1] );
T = nrow(swaps$X );

for( s in 1 : length(Steps) )
{

    # compute series at aggregated time steps
    k = Steps[ s ];
    AggX = NULL;
    t = 1;
    while( ( t + k + 1 ) <= T )
    {
        NewTerm = apply( matrix(swaps$X[ t : (t+k-1), ], ,ncol(swaps$X) ),2,sum);
        AggX = rbind( AggX, NewTerm ); 
        t = t + k;
    }

    # empirical mean/covariance
    
    if(s==1)
    {
        M_hat = matrix(apply(AggX, 2, mean));
        S_hat = cov(AggX);
        Agg[[s]] = list( M_hat = M_hat, S_hat = S_hat, M_norm = k / Steps[ 1 ] * M_hat,  S_norm = k / Steps[ 1 ] * S_hat );
    }else
    {

        Agg[[s]] = list(
                M_hat = matrix(apply(AggX,2,mean)),
                S_hat = cov(AggX),
                # mean/covariance implied by propagation law of risk for invariants
                M_norm = k / Steps[ 1 ] * Agg[[ 1 ]]$M_hat, 
                S_norm = k / Steps[ 1 ] * Agg[[ 1 ]]$S_hat
                );
    }

    # plots
    h1 = TwoDimEllipsoid( Agg[[ s ]]$M_norm, Agg[[ s ]]$S_norm, 1, 0, 0 );
    
    h2 = TwoDimEllipsoid( Agg[[ s ]]$M_hat, Agg[[ s ]]$S_hat, 1, 0, 0 );
}
