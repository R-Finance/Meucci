#' @title Maps moments of log-returns to linear returns .
#'
#' @description Map moments of log-returns to linear returns, as described in  A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   Mu     [vector] (N x 1)
#'	@param   Sigma  [matrix] (N x N)
#'  
#'	@return  M      [vector] (N x 1)
#'  @return  S      [matrix] (N x N)
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#'
#' See Meucci's script for "Log2Lin.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

Log2Lin = function( Mu, Sigma )
{
	M = exp( Mu + (1/2) * diag( Sigma )) - 1;
	S = exp( Mu + (1/2) * diag( Sigma )) %*% t( exp( Mu + ( 1/2 ) * diag(Sigma) ) ) * ( exp( Sigma ) - 1 );

	return( list( M = M, S = S ) );
}