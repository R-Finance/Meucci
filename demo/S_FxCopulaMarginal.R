#' This script displays the empirical copula of a set of market variables, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_FxCopulaMarginal.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export 

### Load data and select the pair to display

library(pracma)
load( "../data/fX.rda" )

Display = c( 1, 2 );  # 1 = Spot USD/EUR; 2 = Spot USD/GBP; 3 = Spot USD/JPY; 

#############################################################################################################
### Define variables (NB: first column is time)

X = apply( log( db_FX$Data[ , 2 : ncol( db_FX$Data ) ] ), 2, FUN = "diff" );

#############################################################################################################
### Compute empirical copula by sorting
NumObs = nrow( X );
K = ncol ( X );

# Sort and get the permutation indices
C = apply( X, 2, "order" ); 


Copula = matrix( NaN, NumObs, K );

for ( k in 1: K)
{
	# scatter plot
    x  = C[ , k ];

    y  = 1 : NumObs;
    xi = 1 : NumObs;
    yi = interp1(x, y, xi);
    Copula[ , k ] = yi / ( NumObs + 1 );
}

############################################################################################################
### Display

# marginals
NumBins = round(10 * log(NumObs));

layout( matrix(c(1,2,2,1,2,2,0,3,3), 3, 3, byrow = TRUE), heights=c(1,2,1));



#hist( X[ , Display[ 2 ] ], NumBins, xlab = db_FX$Fields[[ Display[ 2 ] + 1 ]], ylab = "", main = "");
barplot( table( cut( X[ , Display[ 2 ] ], NumBins )), horiz=TRUE, yaxt="n")
axis( 2, at = seq(0, 100, 20), labels = seq( 0, 1, 0.2 ) );


# scatter plot
plot( Copula[ , Display[ 1 ] ], Copula[ , Display[ 2 ] ], main = "Copula", 
	xlab = db_FX$Fields[[ Display[ 2 ] + 1 ]], ylab = db_FX$Fields[[ Display[ 1 ] + 1 ]] );

hist( X[ , Display[ 1 ] ], NumBins,xlab = db_FX$Fields[[ Display[ 1 ] + 1 ]], ylab = "", main = "");

