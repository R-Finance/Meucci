#' This script computes the mean-variance frontier of a set of options
#' Described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 6.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_MeanVarianceCalls.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}
##################################################################################################################
### Load dat

data("db" );

##################################################################################################################
### Inputs

# market
Stock_0 = db$Stock[ nrow(db$Stock), ];
Vol_0   = db$Vol[ nrow(db$Stock), ];
Strike  = Stock_0; # ATM strike
Est = apply( diff( db$Dates ), 2, mean ) / 252; # estimation interval
Hor = 2 * Est; # investment horizon

# constraints
N = length( Vol_0 );
Constr = list( Aeq = matrix( 1, 1, N ), beq = 1, Aleq = rbind( diag( 1, N ), -diag( 1, N ) ), 
			   bleq = rbind( matrix( 1, N, 1 ), matrix( 0, N, 1 ) ) 
			  );

J = 10000; # num simulations

##################################################################################################################
### Allocation process
# quest for invariance
x_Stock = diff( log( db$Stock ) );
x_Vol   = diff( log( db$Vol ) );

# estimation
M = matrix( apply(cbind( x_Stock, x_Vol ), 2, mean ) );
S = cov( cbind( x_Stock, x_Vol ) );

# projection
M_Hor = M * Hor / Est;
S_Hor = S * Hor / Est;
X = rmvnorm( J, M_Hor, S_Hor, method = "svd" );
X_Stock = X[ , 1:N ];
X_Vol = X[ , (N + 1):ncol(X) ];

Stock_Hor = repmat( Stock_0,  J, 1  ) * exp( X_Stock );
Vol_Hor   = repmat( Vol_0,  J, 1  ) * exp( X_Vol );

##################################################################################################################
### Pricing
Call_0 = NULL;
Call_Hor = NULL;
for( n in 1 : N )
{
    Rate = 0.04;
    Call_0   = cbind( Call_0, BlackScholesCallPrice( Stock_0[ n ], Strike[ n ], Rate, Vol_0[ n ], db$Expiry[ n ] )$c );
    Call_Hor = cbind( Call_Hor, BlackScholesCallPrice(Stock_Hor[ , n ], Strike[ n ], Rate, Vol_Hor[  ,n ], db$Expiry[ n] - Hor )$c );
}

##################################################################################################################
### Mean-variance
L = Call_Hor / repmat( Call_0, J, 1  ) - 1;
ExpectedValues = matrix( apply( L, 2, mean) );
Covariance = cov( L );
NumPortf = 40;
#[e, vol, w] = 
EFR = EfficientFrontierReturns( NumPortf, Covariance, ExpectedValues, Constr );

##################################################################################################################
### Plots
PlotVolVsCompositionEfficientFrontier( EFR$Composition, EFR$Volatility );

