#' Perform the horizon projection of a Student t invariant, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'  @param	nu    : [scalar] degree of freedom
#'	@param	s     : [scalar] scatter parameter
#'  @param	m     : [scalar] location parameter
#'  @param	T     : [scalar] multiple of the estimation period to the invesment horizon 
#'  
#'  @return	x_Hor : [scalar]
#'  @return	f_Hor : [scalar]
#'	@return	F_Hor : [scalar]
#'
#' @references
#' \url{http://}
#' See Meucci's script for "ProjectionStudentT.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

ProjectionStudentT = function(nu, m, s, T)
{
	# set up grid
	N  = 2 ^ 14; # coarseness level
	a  = -qnorm( 10^(-15), 0, sqrt( T ) );
	h  = 2 * a / N;
	Xi = seq(-a+h, a, h );

	# discretized initial pdf (standardized)
	f = 1 / h * ( pt( Xi + h/2, nu ) - pt( Xi - h/2, nu ) );
	f[ N ] = 1 / h *( pt(-a + h/2, nu ) - pt( -a, nu ) + pt( a, nu )- pt( a-h/2, nu ) );

	# discretized characteristic function
	Phi = fft(f);                     

	# projection of discretized characteristic function
	Signs = ( -1 )^((0:(N-1)) * ( T - 1));   
	Phi_T = h ^ ( T - 1 ) * Signs * (Phi ^ T);

	# horizon discretized pdf (standardized)
	f_T = as.numeric( ifft( Phi_T ) );

	# horizon discretized pdf and cdf (non-standardized)
	x_Hor = m * T + s * Xi;
	f_Hor = f_T / s;
	F_Hor = h * cumsum( f_Hor * s );

	return( list( x = x_Hor, f = f_Hor, F = F_Hor ) );

}