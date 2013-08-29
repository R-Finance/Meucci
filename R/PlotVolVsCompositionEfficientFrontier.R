#' Plot the efficient frontier in the plane of portfolio weights versus standard deviation,
#' as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   Portfolios: [matrix] (M x N) of portfolios weights
#'	@param   vol       : [vector] (M x 1) of volatilities
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "PlotVolVsCompositionEfficientFrontier.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

PlotVolVsCompositionEfficientFrontier = function( Portfolios, vol )
{

	colors = c( "cyan","white","magenta","green","black","red" );
	numcolors = length(colors);

	dev.new();
	xx = dim( Portfolios )[ 1 ];
	N  = dim( Portfolios )[ 2 ];
	
	Data = t( apply( Portfolios, 1, cumsum ) );
	plot(c(0,0), xlim= c( min(vol)*100, max(vol)*100), ylim = c(0, max(Data)), xlab = "Risk %", ylab = "Portfolio weights") ;

	for( n in 1 : N )
	{
	    x = rbind( 1, matrix( 1 : xx ), xx );
	    v = rbind( min(vol), vol, max(vol) ) * 100
	    y = rbind( 0, matrix(Data[ , N-n+1 ]), 0 );
	    polygon( v, y, col = colors[ mod( n, numcolors ) + 1 ] );
	}

}