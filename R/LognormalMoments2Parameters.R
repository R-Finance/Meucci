#' @title Computes the mean and standard deviation of a lognormal distribution from its parameters.
#'
#' @description Computes the mean and standard deviation of a lognormal distribution from its parameters, as described in  
#'  A. Meucci, "Risk and Asset Allocation", Springer, 2005.
#'
#' \deqn{\sigma^{2} = \ln \left( 1 + \frac{V}{E^{2}} \right) , }
#' \deqn{\mu = \ln(E) - \frac{1}{2} \ln \left( 1 + \frac{V}{E^{2}} \right) .}
#'
#'
#'	@param  e     [scalar] expected value of the lognormal distribution
#'  @param	v  	  [scalar] variance of the lognormal distribution
#'  
#'  @return	mu 	  [scalar] expected value of the normal distribution
#'  @return	sig2  [scalar] variance of the normal distribution
#'  
#'  @note	Inverts the formulas (1.98)-(1.99) in "Risk and Asset Allocation", Springer, 2005.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, "E 25- Simulation of a lognormal random variable".
#'
#' See Meucci's script for "LognormalMoments2Parameters.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#determines $\mu$ and $\sigma^2$ from $\Expect\{X\}$ and $\Var\{X\}$, and uses it to determine $\mu$ 
#  and $\sigma^{2}$ such that $\Expect\left\{  X\right\} \equiv 3$ and $\Var\left\{  X\right\}  \equiv 5$
LognormalMoments2Parameters = function( e, v )
{
	sig2 = log( 1 + v / ( e^2 ) );
	mu = log( e ) - sig2 / 2;
	
	return( list( sigma_square = sig2 , mu = mu ) );

}