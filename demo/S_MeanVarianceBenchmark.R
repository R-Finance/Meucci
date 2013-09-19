#' This script projects the distribution of the market invariants for the bond and stock markets
#' (i.e. the changes in yield to maturity and compounded returns) from the estimation interval to the investment horizon.
#' Then it computes the distribution of prices at the investment horizon and translates this distribution into the returns 
#' distribution.
#' Finally, it computes the mean-variance efficient frontier both for a total-return and for a benchmark-driven investor
#' Described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 6.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 257 - Benchmark driven allocation I" and "E 258 - Benchmark driven allocation II".
#'
#' See Meucci's script for "S_MeanVarianceBenchmark.m" and "E 255 - Mean-variance pitfalls: two-step approach II" from the book.
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
data("stockSeries"); 	

###################################################################################################################
### Inputs

tau       = 4 / 52; # time to horizon expressed in years
tau_tilde = 1 / 52; # estimation period expressed in years
FlatCurve = 0.04;   
TimeToMat = 5 / 52; # time to maturity of bond expressed in years

# parameters of the distribution of the changes in yield to maturity
u_minus_tau = TimeToMat - tau;
mu = 0 * u_minus_tau;
sigma=( 20 + 5 / 4 * u_minus_tau ) / 10000;

nSim   = 100000;
Budget = 100;

##################################################################################################################
### Estimation of weekly invariants stock market (compounded returns)
Week_C = log( StockSeries$Prices.TimeSeries[ -1, ]) - log( StockSeries$Prices.TimeSeries[-nrow(StockSeries$Prices.TimeSeries), ] );
N = ncol(Week_C);

la = 0.1;
Shrk_Exp = matrix( 0, N, 1 );
Exp_C_Hat = ( 1 - la ) * matrix( apply( Week_C, 2, mean ) ) + la * Shrk_Exp;

lb = 0.1;
Shrk_Cov = diag( 1, N ) * sum(diag(cov(Week_C))) / N;
Cov_C_Hat = ( 1 - lb ) * cov(Week_C) + lb * (Shrk_Cov);

######################################################################################################(############
### Stock market projection to horizon and pricing 
Exp_Hrzn_C_Hat = Exp_C_Hat * tau / tau_tilde;
Cov_Hrzn_C_Hat = Cov_C_Hat * tau / tau_tilde;
StockCompReturns_Scenarios = rmvnorm( nSim, Exp_Hrzn_C_Hat, Cov_Hrzn_C_Hat);

StockCurrent_Prices = matrix( StockSeries$Prices.TimeSeries[ nrow(StockSeries$Prices.TimeSeries), ]);
StockMarket_Scenarios = ( matrix( 1, nSim, 1) %*% t(StockCurrent_Prices)) * exp( StockCompReturns_Scenarios );

##################################################################################################################
# MV inputs - analytical
Stock = ConvertCompoundedReturns2Price(Exp_Hrzn_C_Hat, Cov_Hrzn_C_Hat, StockCurrent_Prices);
print( Stock$Exp_Prices );
print( Stock$Cov_Prices );

##################################################################################################################
# MV inputs - numerical
StockExp_Prices = matrix( apply( StockMarket_Scenarios, 2, mean ));
StockCov_Prices = cov( StockMarket_Scenarios );
print(StockExp_Prices);
print(StockCov_Prices);

##################################################################################################################
### Bond market projection to horizon and pricing 
BondCurrent_Prices_Shifted = exp( -FlatCurve * u_minus_tau);
BondCurrent_Prices = exp( -FlatCurve * TimeToMat);

# generate changes in yield-to-maturity
DY_Scenarios = matrix( rnorm( nSim, mu * tau / tau_tilde, sigma * sqrt(tau / tau_tilde)) ); 
# compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
X = -u_minus_tau * DY_Scenarios;
BondMarket_Scenarios = BondCurrent_Prices_Shifted * exp(X); 

# MV inputs - analytical
Exp_Hrzn_DY_Hat  = mu * tau / tau_tilde;
SDev_Hrzn_DY_Hat = sigma * sqrt(tau / tau_tilde);
Cov_Hrzn_DY_Hat  = diag( SDev_Hrzn_DY_Hat, length(SDev_Hrzn_DY_Hat) ) %*% diag( SDev_Hrzn_DY_Hat, length(SDev_Hrzn_DY_Hat) );
Bond = ConvertChangeInYield2Price(Exp_Hrzn_DY_Hat, Cov_Hrzn_DY_Hat, u_minus_tau, BondCurrent_Prices_Shifted);
print(Bond$Exp_Prices);
print(Bond$Cov_Prices);

# MV inputs - numerical
BondExp_Prices = t( mean( BondMarket_Scenarios ));
BondCov_Prices = cov(BondMarket_Scenarios);

##################################################################################################################
### Put market together and compute returns
Current_Prices   = rbind( StockCurrent_Prices, BondCurrent_Prices);
Prices_Scenarios = cbind( StockMarket_Scenarios, BondMarket_Scenarios )
Rets_Scenarios   = Prices_Scenarios / (matrix( 1, nSim, 1 ) %*% t(Current_Prices)) - 1;
E = matrix(apply( Rets_Scenarios, 2, mean));
S = cov(Rets_Scenarios);

N = ncol(StockSeries$Prices.TimeSeries) + 1;
w_b = matrix( 1, N, 1) / N; # relative benchmar weights

##################################################################################################################
### Portolio optimization
# MV total return quadratic optimization to determine one-parameter frontier of quasi-optimal solutions 
NumPortf = 40;
Ef = EfficientFrontierReturns( NumPortf, S, E );
Rel_ExpectedValue = array( 0, NumPortf );
Rel_Std_Deviation = array( 0, NumPortf );
for( k in 1 : NumPortf )
{
    Rel_ExpectedValue[ k ] = t( Ef$Composition[ k, ] - w_b) %*% E;
    Rel_Std_Deviation[ k ] = sqrt(t( Ef$Composition[ k, ] - w_b ) %*% S %*% ( Ef$Composition[ k,  ] - w_b ) );
}

##################################################################################################################
### Benchmark-relative statistics
# MV benchmark-relative quadratic optimization to determine one-parameter frontier of quasi-optimal solutions
Ef_b = EfficientFrontierReturnsBenchmark(NumPortf, S, E, w_b);
Rel_ExpectedValue_b = array( 0, NumPortf );
Rel_Std_Deviation_b = array( 0, NumPortf );
for( k in 1 : NumPortf )
{
	Rel_ExpectedValue_b[ k ] = t( Ef_b$Composition[ k, ]  - w_b ) %*% E; 
    Rel_Std_Deviation_b[ k ] = sqrt(t( Ef_b$Composition[ k, ] - w_b ) %*% S %*% ( Ef_b$Composition[ k,  ] - w_b ) );
}

##################################################################################################################
### Plots
# frontiers in total return space
dev.new();
plot( Ef$Volatility, Ef$ExpectedValue, type = "l", lwd = 2, col = "blue", xlab = "st.dev. rets.", ylab = "exp.val rets.",
	xlim =c( Ef_b$Volatility[1], Ef_b$Volatility[length(Ef_b$Volatility)] ), ylim = c( min(Ef_b$ExpectedValue), max(Ef_b$ExpectedValue)) );
lines( Ef_b$Volatility , Ef_b$ExpectedValue, type = "l", lwd = 2, col = "red" );
legend( "topleft", 1.9, c( "total ret", "relative" ), col = c( "blue","red" ),
     lty=1, bg = "gray90" );

# frontiers in relative return space
dev.new();
plot( Rel_Std_Deviation, Rel_ExpectedValue, type = "l", lwd = 2, col = "blue", xlab = "TE rets.", ylab = "EOP rets.",
	xlim =c( Rel_Std_Deviation_b[1], Rel_Std_Deviation_b[length(Rel_Std_Deviation_b)] ), ylim = c( min( Rel_ExpectedValue_b ), max( Rel_ExpectedValue_b )) );
lines( Rel_Std_Deviation_b, Rel_ExpectedValue_b, lwd = 2, col = "red" );
legend( "topleft", 1.9, c( "total ret", "relative" ), col = c( "blue","red" ),
     lty=1, bg = "gray90" );
