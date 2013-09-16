#' @title Solve for G that maximises sample r-square of X*G'*B' with X under constraints A*G<=D
#' and Aeq*G=Deq
#'
#' @description Solve for G that maximises sample r-square of X*G'*B' with X under constraints A*G<=D
#' and Aeq*G=Deq (A,D, Aeq,Deq conformable matrices),as described in  A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005.
#'  
#'  @param   X   : [matrix] (T x N)
#'  @param   B   : [matrix] (T x K)
#'  @param   W   : [matrix] (N x N)
#'  @param   A   : [matrix] linear inequality constraints
#'  @param   D   : [matrix]
#'  @param   Aeq : [matrix] linear equality constraints
#'  @param   Deq : [matrix] 
#'  @param   lb  : [vector] lower bound
#'  @param   ub  : [vector] upper bound
#'  
#'  @return   G   : [matrix] (N x K)
#'
#'  @note
#'  Initial code by Tai-Ho Wang 
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' Used in "E 123 – Cross-section factors: generalized cross-section industry factors".
#'
#' See Meucci's script for "MaxRsqCS.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

MaxRsqCS = function(X, B, W, A = NULL, D = NULL, Aeq = NULL, Deq, lb = NULL, ub = NULL)
{
    N = ncol(X);
    K = ncol(B);

    # compute sample estimates
    Sigma_X = (dim(X)[1]-1)/dim(X)[1] * cov(X);


    # restructure for feeding to quadprog 
    Phi = t(W) %*% W;

    # restructure the linear term of the objective function 
    FirstDegree = matrix( Sigma_X %*% Phi %*% B, K * N, );

    # restructure the quadratic term of the objective function
    SecondDegree = Sigma_X;
    
    for( k in 1 : (N - 1) )
    {
        SecondDegree = blkdiag(SecondDegree, Sigma_X);
    }
    
    SecondDegree = t( kron( sqrt(Phi) %*% B, diag( 1, N ) ) ) %*% SecondDegree %*% kron( sqrt( Phi ) %*% B, diag( 1, N ) );

    # restructure the equality constraints 
    if( !length(Aeq)  )
    {
        AEq = Aeq;
    }else
    {
        AEq = blkdiag(Aeq);
        for( k in 2 : K )
        {
            AEq = blkdiag(AEq, Aeq);
        }
    }

    Deq = matrix( Deq, , 1);

    # resturcture the inequality constraints 
    if( length(A) )
    {
        AA = NULL
        for( k in 1 : N )
        {
            AA = cbind( AA, kron(diag( 1, K ), A[ k ] ) ); ##ok<AGROW>
        }
    }else
    {
         AA = A;
    }

    if( length(D))
    {
        D = matrix( D, , 1 );
    }

    # restructure upper and lower bounds
    if( length(lb) )
    {
        lb =  matrix( lb, K * N, 1 );
    }

    if( length(ub) )
    {
        ub = matrix( ub, K * N, 1 );
    }

    # initial guess
    x0 = matrix( 1, K * N, 1 );
    if(length(AA))
    {
        AA = ( AA + t(AA) ) / 2; # robustify
          
    }

    Amat = rbind( AEq, AA);
    bvec = c( Deq, D  );
    
    # solve the constrained generlized r-square problem by quadprog
    #options = optimset('LargeScale', 'off', 'MaxIter', 2000, 'Display', 'none');
    

    b = ipop( c = matrix( FirstDegree ), H = SecondDegree, A = Amat, b = bvec, l = lb , u = ub , r = rep(0, length(bvec)) )
    # reshape for output
    G = t( matrix( attributes(b)$primal, N, ) );

    return( G );
}