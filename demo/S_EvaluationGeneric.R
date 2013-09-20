#' Determine the optimal allocation, as described in A. Meucci "Risk and Asset Allocation", Springer, 2005
#' 
#'  @param  Market          : [struct] market parameters
#'  @param  InvestorProfile : [struct] investor's parameters
#'  
#'  @return Allocation      : [vector] (N x 1)
#'
#' @note
#' 	Compute optimal allocation, only possible if hidden parameters were known: thus it is not a "decision", we call it a "choice"
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 285 - Estimation risk and opportunity cost".
#'
#' See Meucci's script for " EvaluationChoiceOptimal.m" 
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}


EvaluationChoiceOptimal = function( Market, InvestorProfile )
{
	Exp_Prices = diag( Market$CurrentPrices, length(Market$CurrentPrices) ) %*% ( 1 + Market$LinRets_EV );
	Cov_Prices = diag( Market$CurrentPrices, length(Market$CurrentPrices) ) %*% Market$LinRets_Cov %*% diag( Market$CurrentPrices, length(Market$CurrentPrices) );

	S = solve( Cov_Prices ) %*% diag( 1, dim(Cov_Prices) );
	A = (t( Market$CurrentPrices ) %*% S %*% Market$CurrentPrices)[ 1 ]; 
	B = (t( Market$CurrentPrices ) %*% S %*% Exp_Prices)[1]; 

	Gamma = (( InvestorProfile$Budget - InvestorProfile$RiskPropensity * B) / A )[1];
	Allocation = InvestorProfile$RiskPropensity * S %*% Exp_Prices + Gamma[ 1 ] * S %*% Market$CurrentPrices;

	return( Allocation );
}

#' Compute the certainty-equivalent statisfaction index , as described in A. Meucci "Risk and Asset Allocation", 
#' Springer, 2005.
#'
#'  @param  Allocation      : [vector] (N x 1)
#'  @param  Market          : [struct] market parameters
#'  @param  InvestorProfile : [struct] investor s parameters
#'  
#'  @return CertaintyEquivalent : [scalar]
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 285 - Estimation risk and opportunity cost".
#'
#' See Meucci's script for " EvaluationSatisfaction.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}


EvaluationSatisfaction = function( Allocation, Market, InvestorProfile )
{
	CertaintyEquivalent = t(Allocation) %*% diag( Market$CurrentPrices, length( Market$CurrentPrices ) ) %*% (1 + Market$LinRets_EV) - 1 / (2 * InvestorProfile$RiskPropensity) * t( Allocation ) %*% diag( Market$CurrentPrices, length( Market$CurrentPrices )) %*% Market$LinRets_Cov %*% diag( Market$CurrentPrices, length( Market$CurrentPrices )) %*% Allocation ;

    return( CertaintyEquivalent[1] )
}


#' Determine the allocation of the best performer, as described in A. Meucci "Risk and Asset Allocation", 
#' Springer, 2005.
#'
#'  @param  Market          : [struct] market parameters
#'  @param  InvestorProfile : [struct] investors parameters
#'  
#'  @return Allocation      : [vector] (N x 1)
#'
#' @note
#' 	scenario-dependent decision that tries to pick the optimal allocation
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 285 - Estimation risk and opportunity cost".
#'
#' See Meucci's script for "EvaluationDecisionBestPerformer.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

EvaluationDecisionBestPerformer = function( Market, InvestorProfile )
{
	# find index of best performer
	B = which.max( Market$LinRetsSeries[ nrow(Market$LinRetsSeries) , ] ); ##ok<ASGLU>

	# invest in that asset
	I = diag( 1, length(Market$CurrentPrices) );
	Allocation = InvestorProfile$Budget * I[ , B ] / Market$CurrentPrices[ B ];

	return( Allocation );
}


#' Determine the cost of allocation, as described in A. Meucci "Risk and Asset Allocation", Springer, 2005. 
#'
#'  @param  Allocation      : [vector] (N x 1)
#'  @param  Market          : [struct] market parameters
#'  @param 	InvestorProfile : [struct] investor's parameters
#'  
#'  @return C_Plus          : [scalar] cost
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 285 - Estimation risk and opportunity cost".
#'
#' See Meucci's script for "EvaluationDecisionBestPerformer.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

EvaluationCost = function( Allocation, Market, InvestorProfile )
{
	aXi   = t(Allocation) %*% diag( Market$CurrentPrices, length( Market$CurrentPrices ) ) %*% (1 + Market$LinRets_EV);
	aPhia = t(Allocation) %*% diag( Market$CurrentPrices, length( Market$CurrentPrices ) ) %*% Market$LinRets_Cov %*% diag( Market$CurrentPrices, length( Market$CurrentPrices ) ) %*% Allocation;

	C = ( 1 - InvestorProfile$BaR ) * InvestorProfile$Budget - aXi + sqrt(2 %*% aPhia) * erfinv( 2 * InvestorProfile$Confidence - 1);
	C_Plus = max(C, 0);    
	return( C_Plus );

}


#' This script evaluates a generic allocation decision (in this case the "best performer" strategy, which fully  
#' invest the budget in the last period's best performer).
#' It displays the distribution of satisfaction, cost of constraint violation and opportunity cost for each value 
#' of the market stress-test parameters (in this case the correlation).
#' Described in A. Meucci "Risk and Asset Allocation", Springer, 2005, Chapter 8.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 285 - Estimation risk and opportunity cost".
#'
#' See Meucci's script for "S_EvaluationGeneric.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs
NumScenarios       = 1000;  
NumCorrelations    = 5;
Bottom_Correlation = 0;
Top_Correlation    = 0.99;

##################################################################################################################
### Input investor's parameters
InvestorProfile = NULL;
InvestorProfile$Budget         = 10000;
InvestorProfile$RiskPropensity = 30;
InvestorProfile$Confidence     = 0.9;
InvestorProfile$BaR            = 0.2;

##################################################################################################################
### Input market parameters
NumAssets = 10;  
a      = 0.5; # effect of correlation on expected values and volatility (hidden)
Bottom = 0.06; 
Top    = 0.36; 
Step   = (Top - Bottom) / (NumAssets - 1); 
v      = seq( Bottom, Top, Step ) ; # volatility vector
Market = list();
Market$T = 20; # not hidden
Market$CurrentPrices = 10 * array( 1, NumAssets);  # not hidden

##################################################################################################################
Step = (Top_Correlation - Bottom_Correlation) / (NumCorrelations - 1);
Overall_Correlations = seq( Bottom_Correlation, Top_Correlation, Step ); 

Suboptimal = NULL;
Suboptimal$StrsTst_Satisfaction = NULL;
Suboptimal$StrsTst_CostConstraints = NULL;
Suboptimal$StrsTst_OppCost = NULL;
Optimal = NULL;
Optimal$StrsTst_Satisfaction = NULL;

for( t in 1 : length(Overall_Correlations) )
{
    # input the (hidden) market parameters (only correlations, we assume standard deviations and expected values fixed and known)
    Market$St_Devations = ( 1 + a * Overall_Correlations[ t ]) * v;  # hidden
    Market$LinRets_EV   = 0.5 * Market$St_Devations;       # hidden
    
    Correlation = ( 1 - Overall_Correlations[ t ] ) * diag( 1, NumAssets) + Overall_Correlations[ t ] * matrix( 1, NumAssets, NumAssets);
    Market$LinRets_Cov = diag( Market$St_Devations, length(Market$St_Devations) ) %*% Correlation %*% diag( Market$St_Devations, length(Market$St_Devations) )

    
    ##################################################################################################################
    # compute optimal allocation, only possible if hidden parameters were known: thus it is not a "decision", we call it a "choice"
    Allocation           = EvaluationChoiceOptimal( Market, InvestorProfile );
    Satisfaction_Optimal = EvaluationSatisfaction( Allocation, Market, InvestorProfile );
    
    ##################################################################################################################
    # choose allocation based on available information
    StrsTst_TrueSatisfaction = NULL; 
    StrsTst_CostConstraints  = NULL;
    
    for( s in 1 : NumScenarios )
    {        
        # generate scenarios i_T of information I_T
        Market$LinRetsSeries = rmvnorm( Market$T, Market$LinRets_EV, Market$LinRets_Cov ); 
        
        # scenario-dependent decision that tries to pick the optimal allocation
        Allocation       = EvaluationDecisionBestPerformer( Market, InvestorProfile );   
        TrueSatisfaction = EvaluationSatisfaction( Allocation, Market, InvestorProfile );
        CostConstraints  = EvaluationCost( Allocation, Market, InvestorProfile );
        
        StrsTst_TrueSatisfaction = cbind( StrsTst_TrueSatisfaction, TrueSatisfaction ); ##ok<*AGROW>
        StrsTst_CostConstraints  = cbind( StrsTst_CostConstraints, CostConstraints );        
    }
    
    Suboptimal$StrsTst_CostConstraints = rbind( Suboptimal$StrsTst_CostConstraints, StrsTst_CostConstraints );
    Suboptimal$StrsTst_Satisfaction    = rbind( Suboptimal$StrsTst_Satisfaction, StrsTst_TrueSatisfaction );
    Suboptimal$StrsTst_OppCost         = rbind( Suboptimal$StrsTst_OppCost, Satisfaction_Optimal - StrsTst_TrueSatisfaction + StrsTst_CostConstraints );
    Optimal$StrsTst_Satisfaction       = rbind( Optimal$StrsTst_Satisfaction, Satisfaction_Optimal );    
}

##################################################################################################################
### Display
NumVBins = round(10 * log(NumScenarios));

# optimal allocation vs. allocation decision
for( t in 1 : length(Overall_Correlations) )
{
    dev.new(); 
    par( mfrow = c( 3, 1) )
    hist(Suboptimal$StrsTst_Satisfaction[ t, ], NumVBins, main = "satisfaction", xlab ="", ylab = "" );

    hist(Suboptimal$StrsTst_CostConstraints[ t, ], NumVBins, main = "constraint violation cost", xlab ="", ylab = "");

    hist(Suboptimal$StrsTst_OppCost[ t, ], NumVBins, main = "opportunity cost", xlab ="", ylab = "");
}