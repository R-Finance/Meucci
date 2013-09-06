
#' This script models the joint distribution of the yet-to-be realized key rates of the government curve,
#' as described in A. Meucci "Risk and Asset Allocation", Springer, 2005, chapter 3.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_PureResidualBonds.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
load("../data/bondAttribution.rda");


##################################################################################################################
# bondAttribution$B = key rate durations
# bondAttribution$F = key rate weekly changes
# bondAttribution$X = bonds returns net of carry

Dim = dim(bondAttribution$B);

U = 0 * bondAttribution$X;

for( t in 1 : Dim[1] )
{
    U[ t, ] = bondAttribution$X[ t, ] - bondAttribution$F[ t, ] %*% drop(bondAttribution$B[ t, , ]);
}

C = cor(cbind( U, bondAttribution$F ) );


C_U  = C[ 1:Dim[3], 1:Dim[3] ];
C_FU = C[ 1:Dim[3], -(1:Dim[3]) ];

# not systematic-plus-idiosyncratic model
print(C_U);
print(C_FU); 

