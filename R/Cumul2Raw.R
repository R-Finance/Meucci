#' Map cumulative moments into raw moments, as described in A. Meucci "Risk and Asset Allocation",
#' Springer, 2005
#'
#'	@param    ka  : [vector] (length N corresponding to order N) cumulative moments
#'
#'	@return   mu_ : [vector] (length N corresponding to order N) corresponding raw moments
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "Cumul2Raw.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export
Cumul2Raw = function( ka )
{
	N  = length( ka );
	mu_ = ka;

	for( n in 2 : N )
	{
	    #ka[ n ] = mu_[ n ]; Doesn't make sense

	    for( k in 1 : (n-1) )
	    {
	        mu_[ n ] = mu_[ n ] + choose( n-1, k-1 ) * ka[ k ] * mu_[ n-k ]; 
	    }    
	}

	return( mu_ );
}