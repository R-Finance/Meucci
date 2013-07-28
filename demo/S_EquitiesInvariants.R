#' This file performs the quest for invariance in the stock market, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005, chapter 3.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_EquitiesInvariants.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}



##################################################################################################################
### Load daily stock prices from the utility sector in the S&P 500
load("../data/equities.Rda");

##################################################################################################################
### Pick one stock from database
Stock_Index = 20;
P = Equities$Prices[ 632 : nrow( Equities$Prices ), Stock_Index ]; # select data after 1/8 rule

##################################################################################################################
### Quest for invariance
# first invariant
X = P[ -1 ] / P[ -length( P )];
PerformIidAnalysis( 1 : length( X ), X, 'Analysis for X' );

# second invariant
Y = P[ -1 ] / P[ -length( P )];
PerformIidAnalysis(1 : length( Y ), Y, 'Analysis for Y' );

# third invariant
Z = X ^ 2;
PerformIidAnalysis( 1 : length(Z), Z, 'Analysis for Z' );

# fourth invariant
W = P[ 3 : length( P ) ] - 2 * P[ 2: ( length( P ) -1 ) ] + P[ 1 : ( length( P ) -2 ) ];
PerformIidAnalysis( 1 : length( W ), W, 'Analysis for W' );

