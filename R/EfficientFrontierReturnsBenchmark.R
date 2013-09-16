#' Compute the mean-variance efficient frontier (on returns) by quadratic programming, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005
#'
#'  @param  NumPortf       : [scalar] number of portfolio in the efficient frontier
#'  @param  Covariance     : [matrix] (N x N) covariance matrix
#'  @param  ExpectedValues : [vector] (N x 1) expected returns
#'  @param  Benchmark      : [vector] (N x 1) of benchmark weights
#'  @param  Constraints    : [struct] set of constraints. Default: weights sum to one, and no-short positions
#'  
#'  @return ExpectedValue  : [vector] (NumPortf x 1) expected values of the portfolios
#'  @return Volatility     : [vector] (NumPortf x 1) standard deviations of the portfolios
#'  @return Composition    : [matrix] (NumPortf x N) optimal portfolios 
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "EfficientFrontierReturnsBenchmark.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

EfficientFrontierReturnsBenchmark = function(NumPortf, Covariance, ExpectedValues, Benchmark, Constraints = NULL)
{

    NumAssets = ncol(Covariance);

    ##################################################################################################################
    # determine return of minimum-risk portfolio
    FirstDegree  = -Covariance %*% Benchmark;
    SecondDegree = Covariance;
    if( !length(Constraints) )
    {
        Aeq  = matrix( 1, 1, NumAssets);
        beq  = 1;
        A    = -diag( 1, NumAssets);     # no-short constraint
        b    = matrix( 0, NumAssets, 1); # no-short constraint
    }else
    {
        Aeq = Constraints$Aeq;
        beq = Constraints$beq;
        A   = Constraints$Aleq; # no-short constraint
        b   = Constraints$bleq; # no-short constraint
    }    

    Amat = rbind( Aeq, A);
    bvec = rbind( beq, b);

######    MinVol_Weights = quadprog( SecondDegree, -FirstDegree, A, b, Aeq, beq, [], [], x0, options );
    MinVol_Weights = matrix( solve.QP( Dmat = SecondDegree, dvec = -FirstDegree, Amat = -t(Amat), bvec = -bvec, meq = length( beq )  )$solution );
    MinVol_Return  = t( MinVol_Weights ) %*% ExpectedValues;

    ##################################################################################################################
    ### Determine return of maximum-return portfolio
    MaxRet_Return = max(ExpectedValues);
    MaxRet_Index = which( ExpectedValues == max(ExpectedValues) );
    ##################################################################################################################
    ### Slice efficient frontier in NumPortf equally thick horizontal sectors in the upper branch only
    Step = (MaxRet_Return - MinVol_Return) / (NumPortf - 1);
    TargetReturns = seq( MinVol_Return, MaxRet_Return, Step );

    ##################################################################################################################
    ### Compute the NumPortf compositions and risk-return coordinates of the optimal allocations relative to each slice

    # initialization
    Composition   = matrix( NaN, NumPortf, NumAssets);
    Volatility    = matrix( NaN, NumPortf, 1);
    ExpectedValue = matrix( NaN, NumPortf, 1);

    # start with min vol portfolio
    Composition[ 1, ] = t(MinVol_Weights);
    Volatility[ 1 ]     = sqrt(t(MinVol_Weights) %*% Covariance %*% MinVol_Weights);
    ExpectedValue[ 1 ]  = t(MinVol_Weights) %*% ExpectedValues;

    for( i in 2 : (NumPortf - 1) )
    {
        # determine least risky portfolio for given expected return
        AEq = rbind( matrix( 1, 1, NumAssets), t(ExpectedValues) );
        bEq = rbind( 1, TargetReturns[ i ]);
        Amat = rbind( AEq, A);
        bvec = rbind( bEq, b)

        Weights = t( solve.QP( Dmat = SecondDegree, dvec = -FirstDegree, Amat = -t(Amat), bvec = -bvec, meq = length( bEq )  )$solution );
        #Weights = t(quadprog(SecondDegree, FirstDegree, A, b, AEq, bEq, [], [], x0, options));
        
        Composition[ i, ]   = Weights;
        Volatility[ i ]     = sqrt( Weights %*% Covariance %*% t(Weights));
        ExpectedValue[ i ]  = Weights %*% ExpectedValues;
    }

    # add max ret portfolio
    Weights = matrix( 0, 1, NumAssets);
    Weights[ MaxRet_Index ] = 1;
    Composition[ nrow(Composition), ] = Weights;
    Volatility[ length(Volatility) ]   = sqrt(Weights %*% Covariance %*% t(Weights));
    ExpectedValue[ length(ExpectedValue) ]  = Weights %*% ExpectedValues;

    return( list( ExpectedValue = ExpectedValue, Volatility = Volatility, Composition = Composition ) );

}