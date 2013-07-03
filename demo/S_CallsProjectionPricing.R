library(mvtnorm);
library(pracma);

#'This script projects the distribution of the market invariants for the derivatives market
#'Then it computes the distribution of prices at the investment horizon  as described in A. Meucci,
#'"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_CallsProjectionPricing.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}


##################################################################################################################
### Load data

# load 'spot' for underlying and current vol surface, given by
# 'impVol' for different 'days2Maturity' and 'moneyness' (K/S)
load("../data/implVol.Rda");

##################################################################################################################
### Inputs

tau_tilde = 5;    # estimation step (days)
tau = 40;         # time to horizon (days)
Time2Mats = c( 100, 150, 200, 250, 300 );    # current time to maturity of call options in days
Strikes   = c( 850, 880, 910, 940, 970 );    # strikes of call options, same dimension as Time2Mat

r_free = 0.04;    # risk-free rate
J = 10000;       # number of simulations

##################################################################################################################
numCalls   = length( Time2Mats );
timeLength = length( implVol$spot );
numSurfPoints = length( implVol$days2Maturity ) * length( implVol$moneyness );

##################################################################################################################
### Estimate invariant distribution assuming normality
# variables in X are changes in log(spot) and changes in log(imp.vol)
# evaluated at the 'numSurfPoints' points on the vol surface (vectorized).
X = matrix( 0, timeLength - 1, numSurfPoints + 1 );
# log-changes of underlying spot
X[ , 1 ] = diff( log( implVol$spot ) );

# log-changes of implied vol for different maturities
impVolSeries = matrix( implVol$impVol, timeLength, numSurfPoints );
for( i in 1 : numSurfPoints )
{
    X[ , i+1 ] = diff( log( impVolSeries[ , i ] ) );
}

muX = apply( X , 2, mean );
SigmaX = cov( X );

##################################################################################################################
### Project distribution to investment horizon
muX = muX * tau / tau_tilde;
SigmaX = SigmaX * tau / tau_tilde;

##################################################################################################################
### Linearly interpolate the vol surface at the current time to obtain implied vol for the given calls today, and price the calls
spot_T = implVol$spot[ length(implVol$spot) ];
volSurf_T = drop( implVol$impVol[ length(implVol$impVol[, 1, 1] ), , ]);
time2Mat_T = Time2Mats;
moneyness_T = Strikes/spot_T;

impVol_T    = t( InterExtrapolate( volSurf_T,t( rbind( time2Mat_T, moneyness_T )), list( implVol$days2Maturity,implVol$moneyness ) ) );
callPrice_T = BlackScholesCallPrice( spot_T, Strikes, r_free, impVol_T, Time2Mats/252 )$c;

##################################################################################################################
### Generate simulations at horizon
X_ = rmvnorm( J, muX, SigmaX );

# interpolate vol surface at horizon for the given calls
spot_ = spot_T * exp(X_[ , 1 ] );
impVol_ = matrix( 0, J, numCalls);
for( j in 1 : J )
{
    volSurf       = volSurf_T * exp( matrix(X_[ j, -1 ],length( implVol$days2Maturity ),length( implVol$moneyness )));
    time2Mat_     = Time2Mats - tau;
    moneyness_    = Strikes / spot_[ j ];
    impVol_[ j, ] = t( InterExtrapolate( volSurf,t( rbind( time2Mat_T, moneyness_T )), list( implVol$days2Maturity,implVol$moneyness ) ));
}

# price the calls at the horizon
callPrice_ = matrix( 0, J, numCalls );
for( i in 1 : numCalls )
{
    callPrice_[ , i ] = BlackScholesCallPrice( spot_, Strikes[ i ], r_free, impVol_[ , i ], time2Mat_[ i ] / 252 )$c;  
}

m = nrow( callPrice_ );
n = ncol( callPrice_ );
LinearRets = callPrice_ /kronecker( matrix( 1, J, 1), callPrice_T)-1

NumBins = round(10 * log(J));

for( i in 1 : numCalls)
{
    dev.new();
    par( mfrow = c( 2 , 1));
    hist( callPrice_[ , i ], NumBins, xlab = "call price");
    plot(spot_, callPrice_[ ,i ], xlab = "spot price", ylab = "call price" );
}