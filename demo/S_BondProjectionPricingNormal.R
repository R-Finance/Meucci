
#'This script projects the distribution of the market invariants for the bond markets 
#'(i.e. the changes in yield to maturity) from the estimation interval to the investment horizon 
#'Then it computes the distribution of prices at the investment horizon  as described in A. Meucci,
#'"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_BondProjectionPricingNormal.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' 

##################################################################################################################
### Inputs
tau = 1/52;        # time to horizon expressed in years
tau_tilde = 1/52;  # estimation period expressed in years

FlatCurve  = 0.04;   
TimesToMat = c( 1, 5, 10, 52, 520 ) / 52; # time to maturity of selected bonds expressed in years

# parameters of the distribution of the changes in yield to maturity
u_minus_tau = TimesToMat - tau;
mus = 0 * u_minus_tau;
sigmas = ( 20 + 5 / 4 * u_minus_tau ) / 10000;

nSim = 100000;

##################################################################################################################
### Bond market projection to horizon and pricing 
BondCurrent_Prices_Shifted = exp( -FlatCurve * u_minus_tau );
BondCurrent_Prices = exp( -FlatCurve * TimesToMat );

# project bond market to horizon
N = length( TimesToMat ); # number of bonds
U = runif( nSim );
BondMarket_Scenarios = matrix( 0, nSim, N );
for( n in 1 : N )
{
    # generate co-dependent changes in yield-to-maturity
    DY_Scenarios = qnorm( U, mus[ n ] * tau / tau_tilde, sigmas[ n ] * sqrt( tau / tau_tilde ) ); 

    # compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
    X = -u_minus_tau[ n ] * DY_Scenarios;
    BondMarket_Scenarios[ , n ] = BondCurrent_Prices_Shifted[ n ] * exp( X ); 
}

##################################################################################################################
### MV inputs - analytical
Exp_Hrzn_DY_Hat  = mus * tau / tau_tilde;
SDev_Hrzn_DY_Hat = sigmas * sqrt( tau / tau_tilde );
Corr_Hrzn_DY_Hat = matrix( 1, N, N ); # full co-dependence
Cov_Hrzn_DY_Hat  = diag( SDev_Hrzn_DY_Hat ) %*% Corr_Hrzn_DY_Hat %*% diag( SDev_Hrzn_DY_Hat );
Bond = ConvertChangeInYield2Price( Exp_Hrzn_DY_Hat, Cov_Hrzn_DY_Hat, u_minus_tau, BondCurrent_Prices_Shifted );
print( Bond$Exp_Prices );
print( Bond$Cov_Prices );

##################################################################################################################
### MV inputs - numerical
BondExp_Prices = t( apply(BondMarket_Scenarios, 2, mean) );
BondCov_Prices = cov( BondMarket_Scenarios );
print( BondExp_Prices );
print( BondCov_Prices );

### EOF