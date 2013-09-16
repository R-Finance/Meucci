#'  Compute the Black-Scholes price of a European call or put option
#'  as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   spot  : [scalar] spot price of underlying
#'	@param   K     : [scalar] strike of the call optioon
#'	@param   r     : [scalar] risk free rate as a fraction
#'	@param   vol   : [scalar] volatility of the underlying as a fraction
#'	@param   T     : [scalar] time to maturity in years
#'
#'	@return  c     : [scalar] price of European call(s)
#'  @return  p     : [scalar] price of European put(s)
#'	@return  delta : [scalar] delta of the call(s) or put(s)
#'	@return  cash  : [scalar] cash held in a replicating portfolio
#'
#'	@note
#'	Code is vectorized, so the inputs can be vectors or matrices (but sizes must match)
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "BlackScholesCallPrice.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

BlackScholesCallPrice = function( spot, K, r, vol, T )
{
	d1    = ( log( spot / K ) + ( r + vol * vol / 2) * T) / (vol * sqrt(T));
	d2    = d1 - vol * sqrt(T);
	delta = pnorm(d1);
	cash  =  -K * exp( -r * T ) * pnorm( d2 );
	c     = spot * delta + cash;

	return( list( c = c, delta = delta, cash = cash ) );
}

#' @rdname BlackScholesCallPrice
#' @export

BlackScholesPutPrice = function( spot, K, r, vol, T )
{
	d1    = ( log( spot / K ) + ( r + vol * vol / 2) * T) / (vol * sqrt(T));
	d2    = d1 - vol * sqrt(T);
	delta = pnorm( -d1 );
	cash  =  -K * exp( -r * T ) * pnorm( d2 );
	p 	  = -( spot * delta + cash );

	return( list( put = p, delta = delta, cash = cash ) );
}

#' @rdname BlackScholesCallPrice
#' @export

BlackScholesCallPutPrice = function( spot, K, r, vol, T )
{
	d1    = ( log( spot / K ) + ( r + vol * vol / 2) * T) / (vol * sqrt(T));
	d2    =  d1 - vol * sqrt(T);
	cash  =  -K * exp( -r * T ) * pnorm( d2 );
	c =    spot * pnorm(  d1 ) + cash;
	p = -( spot * pnorm( -d1 ) + cash);

	return( list( call = c, put = p, cash = cash ) );
}
