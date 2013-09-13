#' Entropy Pooling Example - Ranking Information script
#'
#' This script performs ranking allocation using the Entropy-Pooling approach by Attilio Meucci, 
#' as it appears in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine, 
#' October 2008, p 100-106.
#'
#' Most recent version of article and MATLAB code available at
#' http://www.symmys.com/node/158
#'
#' @references 
#' A. Meucci, "Fully Flexible Views: Theory and Practice" \url{http://www.symmys.com/node/158}
#' See Meucci script for "RankingInformation/S_MAIN.m"
#' 
#' @author  Xavier Valls \email{flamejat@@gmail.com}

#############################################################################
# Load panel X of joint returns realizations and vector p of respective probabilities
# In real life, these are provided by the estimation process
#############################################################################
data("returnsDistribution");

###########################################################################################################
# compute and plot efficient frontier based on prior market distribution
###########################################################################################################
Options = list();
Options$NumPortf = 20; # number of portfolios in efficient frontier
Options$FrontierSpan = c( 0.3, 0.9 ); # range of normalized exp.vals. spanned by efficient frontier

EF = RIEfficientFrontier( returnsDistribution$X, returnsDistribution$p, Options );
PlotResults( EF$e, EF$Sdev, EF$Composition, EF$Exps );

###########################################################################################################
# process ordering information (this is the core of the Entropy Pooling approach)
###########################################################################################################

# the expected return of each entry of Lower is supposed to be smaller than respective entry in Upper
Lower = 4;  
Upper = 3;
p_ = ViewRanking( returnsDistribution$X, returnsDistribution$p, Lower, Upper )$p_; 

#confidence
c  = 0.5; 
p_ = ( 1 - c ) * returnsDistribution$p + c * p_ ;

###########################################################################################################
# compute and plot efficient frontier based on posterior market distribution
###########################################################################################################

EF_ = RIEfficientFrontier( returnsDistribution$X, p_, Options );
PlotResults( EF_$e, EF_$Sdev, EF_$Composition, EF_$Exps, Lower, Upper );
