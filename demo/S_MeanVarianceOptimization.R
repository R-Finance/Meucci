#' This script projects the distribution of the market invariants for the bond and stock markets 
#' (i.e. the changes in yield to maturity and compounded returns) from the estimation interval to the investment 
#' horizon
#' Then it computes the distribution of prices at the investment horizon and performs the two-step mean-variance 
#' optimization. 
#' Described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 6.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_MeanVarianceHorizon.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
data("stockSeries" );

##################################################################################################################
### Inputs
tau = 4 / 52;        # time to horizon expressed in years
tau_tilde = 1 / 52;  # estimation period expressed in years

FlatCurve  = 0.04;   
TimesToMat = c( 4, 5, 10, 52, 520 ) / 52; # time to maturity of selected bonds expressed in years

# parameters of the distribution of the changes in yield to maturity
u_minus_tau = TimesToMat - tau;
nu  = 8;
mus = 0 * u_minus_tau;
sigmas = ( 20 + 5 / 4 * u_minus_tau ) / 10000;

Num_Scenarios = 100000;
Budget = 100;
Zeta = 10; # risk aversion parameter

##################################################################################################################
### Estimation of weekly invariants stock market (compounded returns)
Week_C = diff( log( StockSeries$Prices.TimeSeries  ) );
T = dim( Week_C )[1];
N = dim( Week_C )[2]

# shrinkage estimator of mean vector
la = 0.1;
Shrk_Exp = matrix( 0, N, 1 );
Exp_C_Hat = (1 - la) * matrix( apply( Week_C, 2, mean) ) + la * Shrk_Exp;

# shrinkage estimator of covariance
lb = 0.1;
Shrk_Cov = diag( 1, N ) * sum( diag( cov( Week_C ) ) ) / N;
Cov_C_Hat = (1 - lb) * cov( Week_C ) + lb * ( Shrk_Cov );

##################################################################################################################
### Stock market projection to horizon and pricing 
Exp_Hrzn_C_Hat = Exp_C_Hat * tau / tau_tilde;
Cov_Hrzn_C_Hat = Cov_C_Hat * tau / tau_tilde;
StockCompReturns_Scenarios = rmvnorm( Num_Scenarios, Exp_Hrzn_C_Hat, Cov_Hrzn_C_Hat );

StockCurrent_Prices   = matrix( StockSeries$Prices.TimeSeries[ ncol(StockSeries$Prices.TimeSeries),  ] );
StockMarket_Scenarios = ( matrix( 1, Num_Scenarios, 1 ) %*% t( StockCurrent_Prices ) ) * exp( StockCompReturns_Scenarios );

##################################################################################################################
### MV inputs
# analytical
CCR2P = ConvertCompoundedReturns2Price( Exp_Hrzn_C_Hat, Cov_Hrzn_C_Hat, StockCurrent_Prices ); 
print( CCR2P$Exp_Prices );
print( CCR2P$Cov_Prices );

# numerical
StockExp_Prices = matrix( apply( StockMarket_Scenarios, 2, mean) );
StockCov_Prices = cov( StockMarket_Scenarios );
print( StockExp_Prices );
print( StockCov_Prices );

##################################################################################################################
### Bond market projection to horizon and pricing 
BondCurrent_Prices_Shifted = exp( -FlatCurve * u_minus_tau );
BondCurrent_Prices = exp( -FlatCurve * TimesToMat );

# project bond market to horizon
N = length( TimesToMat ); # number of bonds

# generate common source of randomness
U = runif( Num_Scenarios);
BondMarket_Scenarios = matrix( 0, Num_Scenarios, N );
for( n in 1 : N )
{
    # generate co-dependent changes in yield-to-maturity
    DY_Scenarios = qnorm( U, mus[ n ] * tau / tau_tilde, sigmas[ n ] * sqrt( tau / tau_tilde ) ); 

    # compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
    X = -u_minus_tau[ n ] * DY_Scenarios;
    BondMarket_Scenarios[ , n ] = BondCurrent_Prices_Shifted[ n ] * exp( X ); 
}

##################################################################################################################
### MV inputs 

# analytical
Exp_Hrzn_DY_Hat  = mus * tau / tau_tilde;
SDev_Hrzn_DY_Hat = sigmas * sqrt(tau / tau_tilde);
Corr_Hrzn_DY_Hat = matrix( 1, N, N ); # full co-dependence
Cov_Hrzn_DY_Hat  = diag(SDev_Hrzn_DY_Hat, length( SDev_Hrzn_DY_Hat)) %*% Corr_Hrzn_DY_Hat %*% diag(SDev_Hrzn_DY_Hat, length( SDev_Hrzn_DY_Hat));
#[BondExp_Prices, BondCov_Prices]
CCY2P = ConvertChangeInYield2Price(Exp_Hrzn_DY_Hat, Cov_Hrzn_DY_Hat, u_minus_tau, BondCurrent_Prices_Shifted);
print( CCY2P$Exp_Prices );
print( CCY2P$Cov_Prices );

# numerical
BondExp_Prices = matrix( apply( BondMarket_Scenarios,2, mean ) );
BondCov_Prices = cov(BondMarket_Scenarios);
print(BondExp_Prices);
print(BondCov_Prices);

##################################################################################################################
### Portolio optimization
# step 1: MV quadratic optimization to determine one-parameter frontier of quasi-optimal solutions ...
E = rbind( StockExp_Prices, BondExp_Prices[ 2 ] );
S = blkdiag( StockCov_Prices, matrix( BondCov_Prices[ 2, 2] ) );
Current_Prices = rbind( StockCurrent_Prices, BondCurrent_Prices[ 2 ] );
Market_Scenarios = cbind( StockMarket_Scenarios, BondMarket_Scenarios[ , 2 ] );

NumPortf = 40;
# frontier with QP (no short-sales)
#[ExpectedValue, EFP$Std_Deviation, EFP$Composition] 

EFP = EfficientFrontierPrices( NumPortf, S, E,Current_Prices, Budget );

# step 2: ...evaluate satisfaction for all EFP$Composition on the frontier ...
Store_Satisfaction = NULL;
for( n in 1 : NumPortf )
{
  Allocation         = matrix( EFP$Composition[ n, ] );
  Objective_Scenario = Market_Scenarios %*% Allocation;
  Utility            = -exp( -1 / Zeta * Objective_Scenario );
  ExpU               = apply( Utility, 2, mean );  
  Satisfaction       = -Zeta * log( -ExpU );
  Store_Satisfaction = cbind( Store_Satisfaction, Satisfaction ); ##ok<AGROW>
}

# ... and pick the best
Optimal_Index = which.max( Store_Satisfaction );
Optimal_Allocation = EFP$Composition[ Optimal_Index, ];
print(Optimal_Allocation);

##################################################################################################################
### Plots
dev.new()

par(mfrow = c( 2, 1 ) );
# rets MV frontier 
h = plot( EFP$Std_Deviation, EFP$ExpectedValue, "l", lwd = 2, xlab = "st.dev. prices", ylab = "exp.val prices",
 xlim = c( EFP$Std_Deviation[1], EFP$Std_Deviation[ length(EFP$Std_Deviation) ]), ylim = c( min(EFP$ExpectedValue), max(EFP$ExpectedValue) ) );


# satisfaction as function of st.deviation on the frontier
h = plot(EFP$Std_Deviation, Store_Satisfaction, "l", lwd = 2, xlab = "st.dev. prices", ylab = "satisfaction",
 xlim = c( EFP$Std_Deviation[1], EFP$Std_Deviation[ length(EFP$Std_Deviation) ]), ylim = c( min(Store_Satisfaction), max(Store_Satisfaction) ) );

