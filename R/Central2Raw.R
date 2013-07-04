#' Map central moments into raw moments
#'
#'	@param   mu  : [vector] (length N corresponding to order N) central moments
#'
#'	@return  mu_ : [vector] (length N corresponding to order N) corresponding raw moments
#'
#' @references
#' \url{http://}
#' See Meucci's script for "Central2Raw.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

Central2Raw = function(mu)
{
	N   = length(mu);
	mu_ = mu;

	for ( n in 2 : N )
	{
	    mu_[ n ] = ( ( -1 ) ^( n+1 ) ) * ( mu[ 1 ] )^(n);
	    for( k in 1 : (n-1) )
	    {
	        mu_[ n ] =  mu_[ n ] + choose( n, k ) * ( (-1) ^ ( n - k + 1 )) * mu_[ k ] * (mu_[ 1 ]) ^ ( n - k);
	    }
	    mu_[ n ] = mu_[ n ] + mu[ n ];
	}

	return( mu_);
}