#' Computes least information kernel smoothing
#' 
#' This script uses Entropy Pooling to compute least information kernel smoothing, as described in 
#' A. Meucci, "Personalized Risk Management: Historical Scenarios with Fully Flexible Probabilities"
#' GARP Risk Professional, Dec 2010, p 47-51
#'
#' @param   Y       Matrix representing the macroeconomic indicator
#' @param   y       scalar reprenting the target to which Y is expected to be close in the Generalized Empirical Distribution
#' @param   h2      N X N matrix 
#'
#' @return  p       list containing the vector of posterior probabilities and information about the optimization performance.
#' 
#' @references 
#' \url{http://www.symmys.com/node/150}
#' See Meucci script for "LeastInfoKernel.m"
#' 
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

LeastInfoKernel = function( Y, y, h2 )
{
    T   = dim(Y)[1];
    N   = dim(Y)[2];
    Aeq = matrix( 1, 1, T );  # constrain probabilities to sum to one...
    beq = 1;
    # ...constrain the first moments...
    Aeq = rbind( Aeq, t(Y) );
    beq = rbind( beq, y );
    
    if( !is.nan(h2) )
    {
        SecMom = h2 + y %*% t( y );  # ...constrain the second moments...
        for( k in 1:N )
        {
            for( l in k:N )
            {
                Aeq = rbind( Aeq, ( Y[ , k ] * Y[ , l ] ) );
                beq = rbind( beq, SecMom[ k, l ] );
            }
        }
    }
    p_0 = matrix( 1, T, 1 ) / T;
    p   = EntropyProg( p_0, matrix(,0,0), matrix(,0,0), Aeq, beq ); # ...compute posterior probabilities
    return( p$p_ );
}

