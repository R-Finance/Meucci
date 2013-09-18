#'This script projects the distribution of the market invariants for the bond markets 
#'(i.e. the changes in yield to maturity) from the estimation interval (Student t assumption)
#'to the investment horizon. Then it computes the distribution of prices at the investment 
#'horizon  as described in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 141 - Fixed-income market: projection of Student t invariants".
#'
#' See Meucci's script for "S_BondProjectionPricingStudentT.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs

tau = 4/52;        # time to horizon expressed in years
tau_tilde = 1/52;  # estimation period expressed in years

FlatCurve  = 0.04;   
TimesToMat = c( 4, 5, 10, 52, 520 ) / 52; # time to maturity of selected bonds expressed in years

# determine the parameters of the distribution of the invariants (changes in yield to maturity)
Periods = tau / tau_tilde; # number of estimation periods until the investment horizon
u_minus_tau = TimesToMat - tau;

nu  = 8;
mus = 0 * u_minus_tau;
sigmas = ( 20 + 5 / 4 * u_minus_tau ) / 10000;
Num_Scenarios = 100000;

##################################################################################################################
### Projection and pricing 
BondCurrent_Prices_Shifted = exp(-FlatCurve * u_minus_tau);
BondCurrent_Prices = exp(-FlatCurve * TimesToMat);

# generate common source of randomness
U = runif( Num_Scenarios );  

N = length( TimesToMat ); # number of bonds
par( mfrow = c( N,1 ));
for( n in 1 : N )
{
    # project bond market to horizon
    Projection = ProjectionStudentT( nu, mus[ n ], sigmas[ n ], Periods);
    
    # generate co-dependent changes in yield-to-maturity
    DY_Scenarios = interp1( Projection$F, Projection$x, U, method = "linear"); 

    # compute the horizon prices, (3.81) in "Risk and Asset Allocation" - Springer
    X = -u_minus_tau[ n ] * DY_Scenarios;
    Z = BondCurrent_Prices_Shifted[ n ] * exp(X); 
    
    # compute and plot linear returns
    L = Z / BondCurrent_Prices[ n ] - 1;  
    
    #for n=1 histogram represents the only bar (not empty)
    hist(L, round(10 * log(Num_Scenarios)), xlab = paste( "Linear returns for bond",  n  ), main = "" );
    
}
