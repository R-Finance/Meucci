#' Plot the marginals of the normal-inverse-Whishart model.
#' Described in A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'  @param	Mu_Simul       : []
#'	@param	InvSigma_Simul : []
#'  @param	Mu_0           : []
#'  @param  T_0            : []
#'  @param  Sigma_0        : []
#'  @param  Nu_0           : []
#'  @param  Legend         : []
#'  
#'  @note Numerically and analytically the marginal pdf of 
#'		- the first entry of the random vector Mu
#'		- the (1,1)-entry of the random matrix inv(Sigma)
#'		when Mu and Sigma are jointly normal-inverse-Wishart: Mu ~ St(Mu_0,Sigma/T_0)
#'                                                            inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "PlotMarginalsNormalInverseWishart.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

PlotMarginalsNormalInverseWishart = function(Mu_Simul, InvSigma_Simul, Mu_0, T_0, Sigma_0, Nu_0, Legend)
{
	NumSimulations = nrow( Mu_Simul );
	NumBins = round( 10 * log( NumSimulations ));

	dev.new();

	#################################################################################################################
	### Mu
	# plot empirical pdf (histogram)
	par( mfrow = c(2,1) );
	h = hist(Mu_Simul[ , 1 ], NumBins, plot= FALSE);
	D = h$mids[ 2 ] - h$mids[ 1 ];
	n= h$counts /( D * NumSimulations );
	plot( h$mids, n, type = "h", main = bquote(paste( .(Legend)," ",  mu)) );
	#barplot( n );

	# superimpose analytical expectation
	points( Mu_0[ 1 ], 0,  pch = 21, bg = "red" );
	

	# superimpose analytical pdf
	x_lo = min(Mu_Simul[ ,1 ]);
	x_hi = max(Mu_Simul[ ,1 ]);
	x_grid = seq( x_lo, x_hi, (x_hi-x_lo)/100 );
	m = Mu_0[ 1 ];
	s = sqrt(Sigma_0[ 1, 1] / T_0 );
	f = 1 / s * dt( (x_grid - m )/s, Nu_0 );
	lines(x_grid, f ,col = "red" );

	#################################################################################################################
	### Sigma
	# plot empirical pdf (histogram)
	h = hist(InvSigma_Simul[  ,1 ], NumBins, plot= FALSE );
	D = h$mids[ 2 ] - h$mids[ 1 ];
	n= h$counts /( D * NumSimulations );
	plot( h$mids, n, type = "h", main = bquote(paste( .(Legend),  " inv(Sigma)")) );
	
	# superimpose analytical expectation
	InvSigma_0=solve(Sigma_0);
	points(InvSigma_0[ 1, 1 ],0, pch = 21, bg = "red" );


	# superimpose analytical pdf
	x_lo = min(InvSigma_Simul[ ,1 ]);
	x_hi = max(InvSigma_Simul[ ,1 ]);
	x_grid = seq( x_lo, x_hi, (x_hi-x_lo)/100 );
	sigma_square = InvSigma_0[ 1, 1] / Nu_0;
	A = Nu_0 / 2;
	B = 2 * sigma_square;
	f = dgamma(x_grid, shape = A, scale = B);
	lines(x_grid, f, col = "red" );
}