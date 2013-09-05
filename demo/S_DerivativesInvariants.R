#' This script performs the quest for invariance in the derivatives market, as described 
#' in A. Meucci,"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_DerivativesInvariants.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export 

##################################################################################################################
### Load implied vol for options on SPX for different time to maturity and moneyness
# Variable name: derivatives
load('../data/derivatives.Rda');

##################################################################################################################
### Simple univariate test
# select the implied vol time series for a specific time to maturity and moneyness 
maturityIndex  = 1;   # 1..6
moneynessIndex = 4;   # 1..7

##################################################################################################################
### Quest for invariance for changes in implied vol and changes in log implied vol

#saving the sequence in a variable for legibility
eachFiveRowsSeq = seq( 1 , length(derivatives$impVol[ , 1, 1 ]), 5 ); 

X = diff( derivatives$impVol[ eachFiveRowsSeq , maturityIndex, moneynessIndex ] );
PerformIidAnalysis( 1:length(X), X, 'Changes in implied vol');

Y = diff(log(derivatives$impVol[ eachFiveRowsSeq , maturityIndex, moneynessIndex ]));
PerformIidAnalysis( 1:length( Y ), Y, 'Changes in log of implied vol' );

##################################################################################################################
### Multivariate test with AR(1) structure

Dim = dim(derivatives$impVol[ eachFiveRowsSeq  , ,  ]);
Z = matrix(log(derivatives$impVol[ eachFiveRowsSeq  , , ] ), Dim[ 1 ], Dim[ 2 ] * Dim[ 3 ]);
# VAR(1) model by least square
X = Z[ -1,  ];
F = cbind(matrix( 1, Dim[ 1 ]-1, 1),  Z[ -length( Z[1, ] ) , ]);
E_XF = t( X ) %*% F / Dim[ 1 ];
E_FF = t( F ) %*% F / Dim[ 1 ];
B = E_XF %*% solve(E_FF);
Eps = X - F %*% t( B ); # residuals

PerformIidAnalysis(1:dim(Eps)[1], Eps[ , 3 ], "VAR(1) residuals");
