#'  @title Fits a multivariate Ornstein - Uhlenbeck process at estimation step tau.
#'
#'  @description Fit a multivariate OU process at estimation step tau, as described in  A. Meucci 
#' "Risk and Asset Allocation", Springer, 2005
#'
#'  @param  Y   : [matrix] (T x N)
#'  @param  tau : [scalar] time step
#'  
#'  @return Mu  : [vector] long-term means
#'  @return Th  : [matrix] whose eigenvalues have positive real part / mean reversion speed
#'  @return Sig : [matrix] Sig = S * S', covariance matrix of Brownian motions
#'
#' @note
#'	o dY_t = -Th * (Y_t - Mu) * dt + S * dB_t where
#'	o dB_t: vector of Brownian motions
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#'
#' See Meucci's script for "FitOrnsteinUhlenbeck.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

FitOrnsteinUhlenbeck = function( Y, tau )
{
	T = nrow(Y);
	N = ncol(Y);

	X    = Y[ -1,  ];
	F    = cbind( matrix( 1, T-1, 1 ), Y[ -nrow(Y), ] );
	E_XF = t(X) %*% F / T;
	E_FF = t(F) %*% F / T;
	B    = E_XF %*% solve( E_FF );
	if( length( B[ , -1 ] ) != 1 )
	{
		Th = -logm( B[ , -1 ] ) / tau;

	}else
	{
		Th = -log( B[ , -1 ] ) / tau;
	}

	Mu = solve( diag( 1, N ) - B[ , -1 ] ) %*%  B[ , 1 ] ;

	U  = F %*% t(B) - X;
	
	Sig_tau = cov(U);

	N = length(Mu);
	TsT = kron( Th, diag( 1, N ) ) + kron( diag( 1, N ), Th );

	VecSig_tau = matrix(Sig_tau, N^2, 1);
	VecSig = ( solve( diag( 1, N^2 ) - expm( -TsT * tau ) ) %*% TsT ) %*% VecSig_tau;
	Sig = matrix( VecSig, N, N );

	return( list( Mu = Mu, Theta = Th, Sigma = Sig ) )
}