#' @title Performs simple invariance (i.i.d.) tests on a time series.
#'
#' @description This function performs simple invariance (i.i.d.) tests on a time series, as described in
#' A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'  @param	Dates : [vector] (T x 1) dates
#'	@param	Data  : [matrix] (T x N) data
#'  @param	Str   : [string]  title for the plot 
#'  
#'  @note it checks the evolution over time
#'
#'   it checks that the variables are identically distributed by looking at the histogram of two subsamples
#'
#'   it checks that the variables are independent by looking at the 1-lag scatter plot
#'
#'   under i.i.d. the location-dispersion ellipsoid should be a circle
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#'
#' See Meucci's script for "PerformIidAnalysis.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

PerformIidAnalysis = function( Dates = dim( Data)[1], Data, Str = "")
{

	##########################################################################################################
	### Time series over time
	dev.new(); 
	plot( Dates, Data, main = Str ); 
	#datetick( 'x', 'mmmyy', 'keeplimits', 'keepticks' );
	

	##########################################################################################################
	### Test "identically distributed hypothesis": split observations into two sub-samples and plot histogram
	Sample_1 = Data[ 1:round(length(Data) / 2) ];
	Sample_2 = Data[(round(length(Data)/2) + 1) : length(Data) ];
	num_bins_1 = round(5 * log(length(Sample_1)));
	num_bins_2 = round(5 * log(length(Sample_2)));
	X_lim = c( min(Data) - .1 * (max(Data) - min(Data)), max(Data) + .1 * (max(Data) - min(Data)));

	dev.new();

	layout( matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = TRUE), heights=c(1,1,1));
	hist(Sample_1, num_bins_1, xlab = "", ylab = "", main = "first half" );
	hist(Sample_2, num_bins_2, xlab = "", ylab = "", main = "second half" );
	
	##########################################################################################################
	### Test "independently distributed hypothesis": scatter plot of observations at lagged times
	

	X = Data[ 1 : length(Data)-1 ];
	Y = Data[ 2 : length(Data) ];
	plot(X, Y, main="changes in implied vol");

	m = cbind( apply( cbind( X, Y ), 2, mean ));
	S = cov( cbind( X, Y ));
	TwoDimEllipsoid( m, S, 2, FALSE, FALSE);

}