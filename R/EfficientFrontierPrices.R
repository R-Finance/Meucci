#' @title Computes the mean-variance efficient frontier (on prices) by quadratic programming
#'
#' @description Compute the mean-variance efficient frontier (on prices) by quadratic programming, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'  @param   NumPortf        [scalar] number of portfolio in the efficient frontier
#'  @param   Covariance      [matrix] (N x N) covariance matrix
#'  @param   ExpectedValues  [vector] (N x 1) expected returns
#'  @param   Current_Prices  [vector] (N x 1) current prices
#'  @param   Budget          [scalar] budget constraint
#'  
#'  @return  ExpectedValue   [vector] (NumPortf x 1) expected values of the portfolios
#'  @return  Std_Deviation   [vector] (NumPortf x 1) standard deviations of the portfolios
#'  @return  Composition     [matrix] (NumPortf x N) optimal portfolios 
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#'
#' See Meucci's script for "EfficientFrontierReturns.m".
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

EfficientFrontierPrices = function( NumPortf, Covariance, ExpectedValues, Current_Prices, Budget )
{

    NumAssets = ncol( Covariance );

	##################################################################################################################
	### Determine exp value of minimum-variance portfolio
	FirstDegree  = matrix( 0, NumAssets, 1 );
	SecondDegree = Covariance;
	Aeq = t( Current_Prices );
	beq = Budget;
	A   = -diag( 1, NumAssets );
	b   = matrix( 0, NumAssets, 1 );
	Amat = rbind( Aeq, A);
    bvec = rbind( beq, b);
	#x0  = Budget / NumAssets * matrix( 1 , NumAssets, 1 );
	MinVol_Allocation = matrix( solve.QP( Dmat = SecondDegree, dvec = -FirstDegree, Amat = -t(Amat), bvec = -bvec, meq = length( beq )  )$solution );
	MinVol_ExpVal = t( MinVol_Allocation ) %*% ExpectedValues;

	##################################################################################################################
	### Determine exp value of maximum-expected value portfolio
	Max_ExpVal = Budget * apply( ExpectedValues / Current_Prices, 2, max );

	##################################################################################################################
	### Slice efficient frontier in NumPortf equally thick horizontal sectors in the upper branch only
	Target_ExpectedValues = MinVol_ExpVal + ( 0 : NumPortf ) * ( Max_ExpVal - MinVol_ExpVal ) / NumPortf;

	##################################################################################################################
	### Compute the NumPortf compositions and risk-return coordinates
	Composition   = matrix( NaN, NumPortf, NumAssets );
	Std_Deviation = matrix( NaN, NumPortf, 1 );
	ExpectedValue = matrix( NaN, NumPortf, 1 );

	Min_ExpectedValue = min(ExpectedValues);
	Max_ExpectedValue = max(ExpectedValues);

	IndexMin = which.min(ExpectedValues);	
	IndexMax = which.max(ExpectedValues);
	
	for( i in 1 : NumPortf )
	{
	    # determine initial condition
	    Matrix = rbind( cbind( Min_ExpectedValue, Max_ExpectedValue ),
	              		cbind( Current_Prices[ IndexMin ], Current_Prices[ IndexMax ] ) );
	          
	    Allocation_0_MinMax = solve( Matrix ) %*% rbind( Target_ExpectedValues[ i ], Budget );
	    
	    Allocation_0 = matrix( 0, NumAssets, 1 );
	    Allocation_0[ IndexMin ] = Allocation_0_MinMax[ 1 ];
	    Allocation_0[ IndexMax ] = Allocation_0_MinMax[ 2 ];
	    
	    # determine least risky portfolio for given expected return
	    AEq = rbind( Aeq, t(ExpectedValues) );
	    bEq = rbind( beq, Target_ExpectedValues[ i ] );
	    Amat = rbind( AEq, A);
        bvec = rbind( bEq, b);

	    #options = optimset('Algorithm', 'medium-scale');
	    Allocation = t( solve.QP( Dmat = SecondDegree, dvec = -FirstDegree, Amat = -t(Amat), bvec = -bvec, meq = length( bEq )  )$solution );
	    
	    # store 
	    Composition[ i, ]   = Allocation;
	    Std_Deviation[ i ]  = sqrt(Allocation %*% Covariance %*% t(Allocation) );
	    ExpectedValue[ i ]  = Target_ExpectedValues[ i ];
	}

	return( list( ExpectedValue = ExpectedValue, Std_Deviation = Std_Deviation, Composition = Composition ) );
}