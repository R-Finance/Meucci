#' Map raw moments into central moments, as described in A. Meucci "Risk and Asset Allocation",
#'  Springer, 2005
#'
#'	@param    mu_ : [vector] (length N corresponding to order N) corresponding raw moments
#'
#'	@return   mu  : [vector] (length N corresponding to order N) central moments
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "Raw2Central.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

Raw2Central = function( mu_ )
{
	N  = length( mu_ );
	mu = mu_;

	for( n in 2 : N )
	{
		mu[ n ] = ( (-1) ^ n ) * ( mu_[ 1 ] )^( n );
	    
	    for( k in 1 : (n-1) )
	    {
	        mu[ n ] = mu[ n ] + choose( n, k ) * ((-1)^(n-k)) * mu_[ k ] * (mu_[ 1 ])^(n-k) ; 
	    }    

	    mu[ n ] = mu[ n ] + mu_[ n ];
	}

	return( mu );
}