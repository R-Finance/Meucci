#' Convert compounded returns to prices for equity-like securities, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'  @param	Exp_Comp_Rets   : [vector] (N x 1) expected values of compounded returns
#'	@param	Cov_Comp_Rets   : [matrix] (N x N) covariance matrix of compounded returns
#'  @param	Starting_Prices : [vector] (N x 1) 
#'  
#'  @return	Exp_Prices    : [vector] (N x 1) expected values of prices
#'  @return	Cov_Prices    : [matrix] (N x N) covariance matrix of prices
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See (6.77)-(6.79) in "Risk and Asset Allocation"-Springer (2005), by A. Meucci
#' See Meucci's script for "ConvertCompoundedReturns2Price.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export


ConvertCompoundedReturns2Price = function(Exp_Comp_Rets, Cov_Comp_Rets, Starting_Prices)
{
	Mu = log(Starting_Prices) + Exp_Comp_Rets;
	Sigma = Cov_Comp_Rets;

	Exp_Prices = exp( Mu + 0.5 * diag( Sigma ) );
	Cov_Prices = exp( Mu + 0.5 * diag( Sigma ) ) %*% t( exp( Mu + 0.5 * diag(Sigma) )) * ( exp( Sigma ) - 1 );

	return( list( Exp_Prices = Exp_Prices, Cov_Prices = Cov_Prices ) );
}