#' Computes the Call Price
#' 
#' Pricing function to apply to each scenario in order to generate the P&L distribution, as described
#' A. Meucci, "Personalized Risk Management: Historical Scenarios with Fully Flexible Probabilities"
#' GARP Risk Professional, Dec 2010, p 47-51
#'
#' @param   P       matrix of prices
#' @param   K       
#' @param   r       risk
#' @param   t       expiry
#' @param   s       volatility
#'
#' @return  C       Prices
#' 
#' @references 
#' \url{http://www.symmys.com/node/150}
#' See Meucci script for "CallPrice.m"
#' 
#' @author Xavier Valls \email{flamejat@@gmail.com}

CallPrice = function( P, K, r, t, s )
{
    d_1 = log( P/K ) + ( r + s * s/2 ) * t;
    d_2 = d_1 - s * sqrt( t );

    C = P * pnorm( d_1 ) - K * exp( -r * t ) * pnorm( d_2 );

    return( C );
}



#'This script uses Entropy Pooling to compute Fully Flexible Probabilities for historical scenarios
#'based on time periods, market conditions, constraints on moments, etc., as described in
#'A. Meucci, "Personalized Risk Management: Historical Scenarios with Fully Flexible Probabilities"
#'GARP Risk Professional, Dec 2010, p 47-51
#'
#' Most recent version of article and code available at
#' http://www.symmys.com/node/150
#' @references 
#' \url{http://www.symmys.com/node/150}
#' See Meucci script for "DoubleDecay.m"
#' 
#' @author Xavier Valls \email{flamejat@@gmail.com}

##########################################################################
# risk drivers scenarios
###########################################################################

load( "../data/dbFFP.rda" )

Infl = dbFFP$Data[ , length( dbFFP$Names ) ];
Vix = dbFFP$Data[ , length( dbFFP$Names ) - 1 ];
Crude = dbFFP$Data[ , length( dbFFP$Names )-3 ];
Swp10 = dbFFP$Data[ , 2 ];
SnP = dbFFP$Data[ , 4 ];

X = diff( log( cbind( SnP, Vix, Swp10 ) ) );
Y = matrix(Infl[ -nrow( dbFFP$Data ) ]);

##########################################################################
#assign probabilities to historical scenarios
###########################################################################
# DefineProbs = "1" : rolling window
# DefineProbs = "2" : exponential smoothing
# DefineProbs = "3" : market conditions
# DefineProbs = "4" : kernel damping
# DefineProbs = "5" : partial information prox. kernel damping
# DefineProbs = "6" : partial information: match covariance

DefineProbs = 1;

T = dim(X)[1];
p = matrix( 0, T, 1 );


if( DefineProbs == 1)
{
	# rolling window

        tau = 2 * 252;
        p[ 1:tau ] = 1;
        p = p / sum( p );
} else if( DefineProbs == 2 )
{ 	
	# exponential smoothing

        lmd = 0.0166;
        p   = exp( -lmd * ( T - ( 1 : T ) ) );
        p   = p / sum( p );

} else if( DefineProbs == 3 )
{ 
	# market conditions
        Cond = Y >= 2.8;
        p[ Cond ] = 1;
        p = p / sum( p );

} else if( DefineProbs == 4 )
{ 
	# kernel damping
        y  = 3;
        h2 = cov( matrix( diff( Y ) ) );
        p  = dmvnorm( Y, y, h2 );
        p  = p / sum( p );
    
} else if( DefineProbs == 5 )
{ 
	# partial information prox. kernel damping
        y  = 3;
        h2 = NaN; # set h2=NaN for no conditioning on second moments
        h2 = cov( 1 * diff( Y ) );
        p  = LeastInfoKernel( Y, y, h2 );
    
} else if( DefineProbs == 6 ){ 
	 #partial information: match covariance

		l_c = 0.0055;
		l_s = 0.0166;

		N = 20;
		Dd = DoubleDecay( X, l_c, l_s );

		p = Fit2Moms( X, Dd$m, Dd$S );
}

###########################################################################
# P&L scenarios
###########################################################################

N = 20;

# call parameters
S_0    = SnP[ length(SnP) ];
vol_0  = Vix[ length(Vix)];
rf_0   = Swp10[ length(Swp10) ];
K      = S_0 * ( seq( 0.8, 1.1, length = N) );
Expiry = ( 2: (N+1) ) / 252;

S_T   = S_0 * exp( X[ , 1 ] );
vol_T = vol_0 * exp( X[ , 2 ] );
rf_T  = rf_0 * exp( X[ , 3 ] );

PnL = matrix( NaN, T, N );

# securities scenarios
for( n in 1:N )
{
    Call_1 = CallPrice( S_T, K[ n ], rf_T, Expiry[ n ] - 1 / 252, vol_T );
    Call_0 = CallPrice( S_0, K[ n ], rf_0, Expiry[ n ], vol_0 );
    PnL[ , n ] = Call_1 - Call_0;
}

# portfolio scenarios
u     = -rbind( -matrix( 1, N/2, 1 ), matrix( 1, N/2, 1 ) ); # number of units (contracts/shares/etc)
PnL_u = PnL %*% u;



