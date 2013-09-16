#' This script compares hedging based on Black-Scholes deltas with Factors on Demand hedging, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005, Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 127 – Factors on demand: no-Greek hedging".
#'
#' See Meucci's script for "S_HedgeOptions.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
data("implVol" );

##################################################################################################################
### Inputs
tau_tilde = 5; # estimation step (days)
tau = 40;      # time to horizon (days)
Time2Mats = cbind( 100, 150, 200, 250, 300 ); # current time to maturity of call options in days
Strikes   = cbind( 850, 880, 910, 940, 970 ); # strikes of call options, same dimension as Time2Mat

r_free = 0.04; # risk-free rate
J = 10000;     # number of simulations

##################################################################################################################
### Underlying and volatility surface
numCalls      = length( Time2Mats );
timeLength    = length( implVol$spot );
numSurfPoints = length( implVol$days2Maturity ) * length( implVol$moneyness );

##################################################################################################################
### Estimate invariant distribution assuming normality
### variables in X are changes in log(spot) and changes in log(imp.vol)
### evaluated at the 'numSurfPoints' points on the vol surface (vectorized).
X = matrix( 0, timeLength-1, numSurfPoints + 1 );
# log-changes of underlying spot
X[ , 1 ] = diff( log( implVol$spot ) );

# log-changes of implied vol for different maturities
impVolSeries = matrix( implVol$impVol, timeLength, numSurfPoints);
for( i in 1 : numSurfPoints )
{
    X[ , i + 1 ] = diff( log( impVolSeries[ , i ] ) );
}
muX = apply(X, 2, mean );
SigmaX = ( dim(X)[1]-1)/dim(X)[1] * cov( X );

##################################################################################################################
### Project distribution to investment horizon
muX = muX * tau / tau_tilde;
SigmaX = SigmaX * tau / tau_tilde;

##################################################################################################################
### Linearly interpolate the vol surface at the current time to obtain
### implied vol for the given calls today, and price the calls
spot_T      = implVol$spot[ length(implVol$spot) ];
volSurf_T   = drop( implVol$impVol[ dim( implVol$impVol )[1], , ] );
time2Mat_T  = Time2Mats;
moneyness_T = Strikes / spot_T;
impVol_T    = t(InterExtrapolate( volSurf_T,cbind( t(time2Mat_T), t( moneyness_T)), list( implVol$days2Maturity, implVol$moneyness ) ) ); # function by John D'Errico
callPrice_T = BlackScholesCallPrice( spot_T, Strikes, r_free, impVol_T, Time2Mats/252 )$c;

##################################################################################################################
### Generate simulations at horizon
X_ = rmvnorm( J, muX, SigmaX );

##################################################################################################################
### Interpolate vol surface at horizon for the given calls
spot_ = spot_T * exp(X_[ , 1 ]);
impVol_ = matrix( 0, J, numCalls);
for( j in 1:J )
{
    volSurf       = volSurf_T *exp( matrix( X_[ j, -1 ], length( implVol$days2Maturity ), length( implVol$moneyness ) ) );
    time2Mat_     = Time2Mats - tau;
    moneyness_    = Strikes / spot_[ j ];
    impVol_[ j, ] = t( InterExtrapolate( volSurf, cbind( t( time2Mat_ ), t( moneyness_) ), list( implVol$days2Maturity, implVol$moneyness ) ) );  # function by John D'Errico
}

##################################################################################################################
### Price the calls
callPrice_ = matrix( 0, J, numCalls );
for( i in 1 : numCalls )
{
    callPrice_[ , i ] = BlackScholesCallPrice( spot_, Strikes[ i ], r_free, impVol_[ , i ], time2Mat_[ i ] / 252 )$c;
}
# linear returns of the calls
Rc = callPrice_ /  repmat( callPrice_T, J, 1) - 1 ;
# linear returns of the underlying
Rsp = spot_ / spot_T - 1;

##################################################################################################################
# Compute the OLS linear (affine) model: Rc = a + b*Rsp + U
Z = cbind( array( 1, J), Rsp );
olsLoadings = ( t(Rc) %*% Z) %*% solve( t(Z) %*% Z );
a = olsLoadings[ , 1 ];
b = olsLoadings[ , 2 ];

##################################################################################################################
# Compute Black-Scholes delta and cash held in replicating portfolio
BSCP = BlackScholesCallPrice( spot_T, Strikes, r_free, impVol_T, Time2Mats / 252 );
a_bs = BSCP$cash / BSCP$c * r_free * tau / 252;
b_bs = t( BSCP$delta / BSCP$c * spot_T);

printf( "OLS: a = [ %s\t]\n", sprintf("\t%7.4f", t(a) ) );
printf( "B-S: a = [ %s\t]\n", sprintf("\t%7.4f", t(a_bs) ) );
printf( "OLS: b = [ %s\t]\n", sprintf("\t%7.4f", t(b) ) );
printf( "B-S: b = [ %s\t]\n", sprintf("\t%7.4f", t(b_bs) ) );

for( i in 1 : numCalls )
{
    dev.new();
    plot( Rsp, Rc[ , i ], xlab = "return underlying" , ylab = "return call option" );
}

