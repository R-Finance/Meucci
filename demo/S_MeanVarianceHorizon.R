#' This script projects the distribution of the market invariants for stock market (i.e. compounded returns) 
#' from the estimation interval to the investment horizon. 
#' Then it computes the distribution of prices at the investment horizon and performs the two-step mean-variance
#' optimization in terms of returns and relative portfolio weights.
#' Described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 6.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 256 – Mean-variance pitfalls: horizon effect".
#'
#' See Meucci's script for "S_MeanVarianceHorizon.m" and "E 255 - Mean-variance pitfalls: two-step approach II" from the book.
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
data("stockSeries"); 	

##################################################################################################################
### Inputs
tau = 1/252;        # time to horizon expressed in years
tau_tilde = 1/52;  # estimation period expressed in years
nSim   = 10000;
Budget = 100;
Zeta   = 10; # risk aversion parameter

##################################################################################################################
### Estimation of weekly invariants stock market (compounded returns)
Week_C = diff( log( StockSeries$Prices.TimeSeries ) );
N = ncol( Week_C );

la = 0.1;
Shrk_Exp = matrix( 0, N, 1);
Exp_C_Hat = (1 - la) * matrix( apply( Week_C, 2, mean ) )  + la * Shrk_Exp;

lb = 0.1;
Shrk_Cov = diag( 1, N ) * sum( diag( cov( Week_C ) ) ) / N;
Cov_C_Hat = (1 - lb) * cov(Week_C) + lb * (Shrk_Cov);

##################################################################################################################
### Stock market projection to horizon and pricing 
Exp_Hrzn_C_Hat = Exp_C_Hat * tau / tau_tilde;
Cov_Hrzn_C_Hat = Cov_C_Hat * tau / tau_tilde;
StockCompReturns_Scenarios = rmvnorm( nSim, Exp_Hrzn_C_Hat, Cov_Hrzn_C_Hat);

StockCurrent_Prices = matrix( StockSeries$Prices.TimeSeries[ nrow( StockSeries$Prices.TimeSeries ), ]);
StockMarket_Scenarios = ( matrix( 1, nSim, 1) %*% t(StockCurrent_Prices)) * exp( StockCompReturns_Scenarios );
##################################################################################################################
### MV inputs - analytical
Stock = ConvertCompoundedReturns2Price(Exp_Hrzn_C_Hat, Cov_Hrzn_C_Hat, StockCurrent_Prices);
print( Stock$Exp_Prices );
print( Stock$Cov_Prices );

##################################################################################################################
### MV inputs - numerical
StockExp_Prices = matrix( apply( StockMarket_Scenarios, 2, mean ));
StockCov_Prices = cov( StockMarket_Scenarios );
print(StockExp_Prices);
print(StockCov_Prices);

StockExp_LinRets = StockExp_Prices / StockCurrent_Prices - 1;
StockCov_LinRets = diag( c(1 / StockCurrent_Prices) ) %*% StockCov_Prices %*% diag( c(1 / StockCurrent_Prices) );

##################################################################################################################
### Portolio optimization
# step 1: MV quadratic optimization to determine one-parameter frontier of quasi-optimal solutions ...
NumPortf = 40;
EFR = EfficientFrontierReturns( NumPortf, StockCov_LinRets, StockExp_LinRets );

# step 2: ...evaluate satisfaction for all allocations on the frontier ...
Store_Satisfaction = NULL;
for( n in 1 : NumPortf )
{
  Allocation         = matrix( EFR$Composition[ n, ] ) * Budget / StockCurrent_Prices;
  Objective_Scenario = StockMarket_Scenarios %*% Allocation;
  Utility            = -exp( -1 / Zeta * Objective_Scenario);
  ExpU               = apply( Utility, 2, mean );  
  Satisfaction       = -Zeta * log( -ExpU );
  Store_Satisfaction = cbind( Store_Satisfaction, Satisfaction ); ##ok<AGROW>
}

# ... and pick the best
Optimal_Index 	   = which.max(Store_Satisfaction);
Optimal_Allocation = EFR$Composition[ Optimal_Index, ];

##################################################################################################################
### Plots
dev.new();

par(mfrow = c( 2, 1 ) );
# rets MV frontier 
h = plot(EFR$Volatility, EFR$ExpectedValue, "l", lwd = 2, xlab = "st.dev. rets.", ylab = "exp.val rets.",
 xlim = c( EFR$Volatility[1], EFR$Volatility[ length(EFR$Volatility) ]), ylim = c( min( EFR$ExpectedValue ), max( EFR$ExpectedValue ) ) );


# satisfaction as function of st.deviation on the frontier
h = plot( EFR$Volatility, Store_Satisfaction, "l", lwd = 2, xlab = "st.dev. rets.", ylab = "satisfaction",
 xlim = c( EFR$Volatility[1], EFR$Volatility[ length(EFR$Volatility) ]), ylim = c( min(Store_Satisfaction), max(Store_Satisfaction) ) );

