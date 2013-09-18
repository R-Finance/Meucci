#' @title Compute recursively the ML estimators of location and scatter of a multivariate Student t distribution with 
#' given degrees of freedom.
#' 
#' @description Compute recursively the ML estimators of location and scatter of a multivariate Student t distribution with 
#' given degrees of freedom, as described in  A. Meucci, "Risk and Asset Allocation", Springer, 2005, section 4.3 - "Maximum likelihood estimators"
#' 
#' Returns \deqn{ \widehat{\mu} = \sum\limits_{t= 1}^{T} \frac{w_{t}}{\sum\limits_{s= 1}^{T}w_{s}} x_{t} ,} 
#' \deqn{\widehat{\Sigma} = \frac{1}{T}\sum\limits_{t= 1}^{T}w_{t}(x_{t}-\widehat{\mu})(x_{t}-\widehat{\mu})^\prime} 
#' 
#' where, adapted to the Sudent T distribution, \deqn{ w_{t}\equiv \frac{\nu+N}{\nu+(x_{t}-\widehat{\mu})^\prime \widehat{\Sigma}^{-1}(  x-\widehat{\mu})  }}
#'
#'  @param   x         : [matrix] (T x N) observations
#'  @param   Nu        : [scalar] degrees of freedom parameter
#'  @param   Tolerance : [scalar] tolerance parameter. Default: 10^(-10)
#'
#'  @return  Mu        : [vector] (N x 1) mean
#'  @return  Sigma     : [matrix] (N x N) covariance
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 188 - Maximum likelihood estimation of a multivariate Student t distribution".
#'
#' See Meucci's script for "MleRecursionForStudentT.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

MleRecursionForStudentT = function(x, Nu, Tolerance = 10^(-10) )
{

    T      = nrow( x );
    N      = ncol( x );
    Ones_N = matrix( 1, 1, N ); # fixed for fast matrix operation
    Ones_T = matrix( 1, T, 1 ); # fixed for fast matrix operation 

    # initialize variables
    w     = matrix( 1, T, 1 );
    Mu    = matrix( 0, N, 1 );
    Sigma = matrix( 0, N, N );

    Error = 10^6;
    # start main loop
    while( Error > Tolerance )
    {
        
        # update
        Mu_Old    = Mu;
        Sigma_Old = Sigma;
        
        # Step 1
        W  = w %*% Ones_N;
        Mu = matrix( apply(  W * x, 2, sum  ) ) / sum( w );
        
        x_c   = x - Ones_T %*% t(Mu);
        Sigma = t( W * x_c ) %*% x_c / T;

        # Step 2
        InvS = solve(Sigma);
        Ma2  = apply( ( x_c %*% InvS ) * x_c, 1, sum );
        w    = ( Nu + N) / ( Nu + Ma2 );

        # convergence
        Error = sum( diag( (Sigma - Sigma_Old) ^2) ) / N + t(Mu - Mu_Old) %*% ( Mu - Mu_Old ) / N;    
    }

    return( list( Mu = Mu, Sigma = Sigma)  );
}