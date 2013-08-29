#' Plot the efficient frontier, as described in  A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   Portfolios : [matrix] (M x N) M portfolios of size N (weights)
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "PlotCompositionEfficientFrontier.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

PlotCompositionEfficientFrontier = function(Portfolios)
{
	dev.new();

	xx = dim( Portfolios )[ 1 ];
	N  = dim( Portfolios )[ 2 ];
	Data = t( apply( Portfolios, 1, cumsum ) );

	plot( c(2000, 2000), xlim= c( 1, xx ), ylim = c( 0, max(Data) ), xlab = " Portfolio # risk propensity", ylab = "Portfolio composition" );
	
	for( n in 1 : N )
	{
	    x = rbind( 1, matrix(1 : xx), xx );
	    y = rbind( 0, matrix( Data[ , N-n+1 ] ), 0 );
	    polygon( x, y, col = rgb( 0.9 - mod(n,3)*0.2, 0.9 - mod(n,3)*0.2, 0.9 - mod(n,3)*0.2) );
	}

}