#' Compute the Black-Scholes price of a European call option
#' as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   spot  : [scalar] spot price of underlying
#'	@param   K     : [scalar] strike of the call optioon
#'	@param   r     : [scalar] risk free rate as a fraction
#'	@param   vol   : [scalar] volatility of the underlying as a fraction
#'	@param   T     : [scalar] time to maturity in years
#'
#'	@return  c     : [scalar] price of European call(s)
#'	@return  delta : [scalar] delta of the call(s)
#'	@return  cash  : [scalar] cash held in a replicating portfolio
#'
#'	@note
#'	Code is vectorized, so the inputs can be vectors or matrices (but sizes must match)
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "BlackScholesCallPrice.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

BlackScholesCallPrice = function(spot, K, r, vol, T)
{
	d1 = ( log( spot / K ) + ( r + vol * vol / 2) * T) / (vol * sqrt(T));
	d2 = d1 - vol * sqrt(T);
	delta = pnorm(d1);
	cash =  -K * exp( -r * T ) * pnorm( d2 );
	c = spot * delta + cash;

	return( list( c = c, delta = delta, cash = cash ) );
}