#' This file performs the quest for invariance in the fixed income market, as described in A. Meucci 
#' "Risk and Asset Allocation", Springer, 2005, Chapter 3.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_FixedIncomeInvariants.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load government yield curve and bond yield data for different dates
data("fixedIncome");

##################################################################################################################
### Pick time-to-maturity for one point on the yield curve
ycMaturityIndex = 4; # 1..6

# select the yield time series for a constant time-to-maturity
yield = fixedIncome$ycYieldPercent[ , ycMaturityIndex ]; 

##################################################################################################################
### Quest for invariance
# changes in the yield curve
X = yield[ -1 ] - yield[ -length( yield ) ];
PerformIidAnalysis( 1:length( X ), X, "Changes in yield curve" );

# changes in the logarithm of the yield curve
Y = log( yield[ -1 ] ) - log( yield[ -length( yield ) ] );
PerformIidAnalysis( 1 : length( Y ), Y, "Changes in log of yield curve" );
