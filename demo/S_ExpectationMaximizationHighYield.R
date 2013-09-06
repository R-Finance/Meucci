#' This script implements the Expectation-Maximization (EM) algoritm, which estimates the parameters of a 
#' multivariate normal distribution when some observations are randomly missing, as described in A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_ExpectationMaximizationHighYield.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
load("../data/highYieldIndices.rda");

##################################################################################################################
### Compute invariants and set NaN for large values
N  = ncol(highYieldIndices$Data);
Series = log( highYieldIndices$Data[ -1, ] ) - log( highYieldIndices$Data[ -nrow(highYieldIndices$Data), ] );
NANs_Index = which( abs( Series ) > (10 ^ 10) );
Series[ NANs_Index ] = NaN;

##################################################################################################################
### Run EM algorithm
FEM = FitExpectationMaximization(Series);

##################################################################################################################
### Display results
dev.new();
par( mfrow = c( N, 1 ) )
for( n in 1 : N )
{
    Drop = is.nan( Series[  , n ] );
    Bad_Dates = highYieldIndices$Dates[ Drop ];
    
    Keep = !is.nan( Series[ , n ] );
    Good_Dates = highYieldIndices$Dates[ Keep ];
    
    plot(Good_Dates[1:(length(Good_Dates)-1)], Series[ Keep, n ], xaxt = "n", xlab = "", ylab = "");
    points(Bad_Dates, FEM$Recovered_Series[ Drop, n ], pch = 21, bg = "red");
    #axis( 1, at = highYieldIndices$Dates, labels=format(highYieldIndices$Dates,"%m/%d/%y"))  # Format x-axis
    
}
legend( "bottom", 1.9, "EM-recovered data",  pch = 21, col = "red" ,bg = "gray90" );


