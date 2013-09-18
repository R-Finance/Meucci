#' This script projects the distribution of the market invariants for the stock market (i.e. the compounded returns)
#' from the estimation interval (normal assumption) to the investment horizon. Then it computes the distribution of prices
#' at the investment horizon analytically, by full Monte Carlo, and by delta/duration approximation.
#'
#' Described in A. Meucci "Risk and Asset Allocation", Springer, 2005, chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 138 - Equity market: linear vs. compounded returns projection II".
#'
#' See Meucci's script for "S_EquityProjectionPricing.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

#################################################################################################################
### Inputs
tau_tilde = 1 / 52;  # estimation period expressed in years
sig = 0.4;
P_T = 1;
NumScenarios = 1000000;
taus    = c( 1/252, 1/52, 1/12, 1, 2 );    # times to horizon expressed in years
tauName = c( '1 day','1 week', '1 month', '1 year', '2 years');

#################################################################################################################
### loop projection/pricing over different times to horizon
for( k in 1 : length( taus ) )
{ 
    tau = taus[ k ];

    # exact simulation of horizon prices
    C_Ttau = rnorm( NumScenarios, 0, sqrt( sig * sig * tau));
    P_Ttau = P_T * exp( C_Ttau );

    # compute analytical pdf
    p_lo = min( P_Ttau );
    p_hi = max( P_Ttau );
    p = seq( p_lo, p_hi, ( p_hi - p_lo ) / 1000 );
    m = log( P_T );
    s = sqrt( sig * sig * tau );
    f = dlnorm( p, m, s );

    # compute approximate pdf
    f_approx = dnorm( p, P_T, sqrt(P_T * P_T * sig * sig * tau));

    # plots
    dev.new()
    NumBins = round(10 * log( NumScenarios));
    hist( P_Ttau, NumBins, freq = FALSE, xlab = "price at the horizon", ylab = "pdf", col="blue",
     main = expression(paste("time to horizon ", tau, " = ", tauName[k])));
    lines( p, f, col = "red");
    lines( p, f_approx, col = "green");
    legend( "topright", 1.9, c( "full Monte Carlo", "analytical", "delta/duration" ), col = c( "blue","red", "green" ),
     lty=1,lwd =c( 10,1,1), bg = "gray90" );
}