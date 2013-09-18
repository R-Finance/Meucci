
#' This script projects summary statistics to arbitrary horizons, as described in A. Meucci 
#' "Risk and Asset Allocation", Springer, 2005, chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 97 - Projection of skewness, kurtosis, and all standardized summary statistics".
#'
#' See Meucci's script for "S_ProjectSummaryStatistics.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs

N = 6;   # focus on first N standardized summary statistics
K = 100; # projection horizon

# generate arbitrary distribution
J = 100000;  # number of scenarios

Z = rnorm( J ); 
X = sin( Z ) + log( cos( Z ) + 2 );

##################################################################################################################
### Compute single-period standardized statistics and central moments
CaSS = CentralAndStandardizedStatistics( X, N );
print( CaSS$ga );
print( CaSS$mu );

# compute single-period non-central moments
mu_ = Central2Raw( CaSS$mu );
print( mu_);

# compute single-period cumulants
ka = Raw2Cumul(mu_);
print(ka);

# compute multi-period cumulants
Ka = K * ka;
print(Ka);

# compute multi-period non-central moments
Mu_ = Cumul2Raw(Ka);
print(Mu_);

# compute multi-period central moments
Mu = Raw2Central(Mu_);
print(Mu);

# compute multi-period standardized statistics
Ga = Mu;
Ga[ 2 ] = sqrt( Mu[ 2 ]);

for( n in 3 : N )
{
    Ga[ n ] = Mu[ n ] / ( Ga[ 2 ] ^ n );
}

print(Ga);