#' @title Plot prior and posterior distributions.
#'
#' @description Plot prior and posterior distributions, as described in  A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005.
#'  
#'  @param   X  : [matrix] (J x N) simulations
#'  @param   p  : [vector] (J x 1) prior probabilities
#'  @param   p_ : [vector] (J x 1) posterior probabilities
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#'
#' See Meucci's script for "pHistPriorPosterior.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

pHistPriorPosterior = function( X, p, p_)
{
    X = as.matrix(X);
    J = dim(X)[1];
    N = dim(X)[2];
    NBins  = round(10 * log(J));

    for( n in 1 : N )
    {
        dev.new();
        
        # set ranges
        xl = min(X[ , n]);
        xh = max(X[ , n]);
        
        par( mfrow = c( 2, 1 ) );

        # prior numerical
        pHist( X[ , n ], p, NBins);
        # xlim([xl, xh]);
        # y1 = ylim();
        # title('prior');
        
        # posterior numerical
        pHist( X[ , n ], p_, NBins);
        # xlim([xl, xh]);
        # y2 = ylim();
        # ylim([min(y1(1), y2(1)), max(y1(2), y2(2))]);
        # title('posterior');
        
        # subplot(2, 1, 1);
        # ylim([min(y1(1), y2(1)), max(y1(2), y2(2))]);
    }

}