#' @title Expectation-Maximization (EM) algorithm to recover missing observations in a time series.
#'
#' @description Expectation-Maximization (EM) algorithm to recover missing observations in a time series ,
#'  as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005, section 4.6.2 "Missing data".
#'  
#'  @param   X         : [matrix] (T x N) of data
#'  
#'  @return  E_EM      : [vector] (N x 1) expectation
#'  @return  S_EM      : [matrix] (N x N) covariance matrix
#'  @return  Y         : [matrix] (T x N) updated data
#'  @return  CountLoop : [scalar] number of iterations of the algorithm
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 177 - Expectation-Maximization algorithm for missing data: formulas"
#' See Meucci's script for "FitExpectationMaximization.m"
#'
#' Dempster, A. P. and Laird, M. N. and Rubin, D. B. - "Maximum Likelihood from Incomplete Data Via the EM Algorithm", 
#' Journal of the Royal Statistical Society, 1977 vol 39 pag. 1-22.
#'
#' Bilmes, J. A.- "A Gentle Tutorial of the EM Algorithm and its Application to Parameter Estimation for Gaussian Mixture
#' and Hidden Markov Models", 1998.
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

FitExpectationMaximization = function(X)
{
    T = nrow(X);
    N = ncol(X);

    # E-M initialization
    idx = apply( !is.nan( X ), 1, all );
    X_Init = X[ idx, ];

    M =  matrix(apply( X_Init, 2, mean ));

    S = cov( X_Init );

    Tolerance = 10 ^ ( -6 ) * mean(rbind( M, sqrt( matrix(diag( S ) ) ) ) );

    # E-M loop
    Convergence = 0;
    CountLoop   = 0;
    Y = X;
    while( !Convergence )
    {
        CountLoop = CountLoop + 1;

        # Step 1: estimation
        C = array( 0, dim = c( T, N, N) );
        for( t in 1 : T )
        {
            Miss = is.nan( X[ t,  ] );
            Obs  = !Miss;
            c = matrix(0, N, N );
            y =  matrix( X[ t,  ]);
            if( any( Miss ) )
            {
                y[ Miss ] = M[ Miss ] + S[ Miss, Obs] %*% ( solve(S[ Obs, Obs ]) %*% matrix (y[ Obs ] - M[ Obs ]) );
                c[ Miss, Miss ] = S[ Miss, Miss ] - S[ Miss, Obs ] %*% ( solve(S[ Obs, Obs ]) %*% S[ Obs, Miss ] );          
            }
    	    Y[ t, ] = y;
            C[ t, , ] = c + (y - M) %*% t(y - M);
        }

        # Step 2: update
        M_new = matrix( apply( Y, 2, mean ));
        S_new = drop( apply( C, c( 2, 3 ), mean ) );
        
        D4 = rbind( ( M_new - M ) ^ 4, matrix(diag( (S_new - S) ^ 2 ) ) );
        Distance = mean( D4 ^ (1/4) );
        Convergence = ( Distance < Tolerance );
        
        M = M_new;
        S = S_new;
    }

    E_EM = M;
    S_EM = S;

    return( list( E_EM = E_EM, S_EM = S_EM, Recovered_Series = Y, CountLoop = CountLoop ) );
}