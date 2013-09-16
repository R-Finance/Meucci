#' Solve for B that maximises sample r-square of F'*B' with X under constraints A*G<=D
#' and Aeq*G=Deq (A,D, Aeq,Deq conformable matrices),as described in  A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005.
#'  
#'  @param   X   : [matrix] (T x N)
#'  @param   F   : [matrix] (T x K)
#'  @param   W   : [matrix] (N x N)
#'  @param   A   : [matrix] linear inequality constraints
#'  @param   D   : [matrix]
#'  @param   Aeq : [matrix] linear equality constraints
#'  @param   Deq : [matrix] 
#'  @param   lb  : [vector] lower bound
#'  @param   ub  : [vector] upper bound
#'  
#'  @return   B   : [matrix] (N x K)
#'
#'  @note
#'  Initial code by Tai-Ho Wang 
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "MaxRsqTS.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

MaxRsqTS = function(X, F, W, A = NULL, D = NULL, Aeq = NULL, Deq, lb = NULL, ub = NULL)
{

	N = dim(X)[ 2 ];
	K = dim(F)[ 2 ];
    X = matrix(as.numeric(X), dim(X))


	# compute sample estimates
	# E_X = apply( X, 2, mean);
	# E_F = apply( F, 2, mean);
	XF = cbind( X, F );
	SigmaJoint_XF = (dim(XF)[1]-1)/dim(XF)[1] * cov(XF);
	
	Sigma_X  = SigmaJoint_XF[ 1:N, 1:N ];
	Sigma_XF = SigmaJoint_XF[ 1:N, (N+1):(N+K) ];
	Sigma_F  = SigmaJoint_XF[ (N+1):(N+K), (N+1):(N+K) ];


    # restructure for feeding to quadprog 
    Phi = t(W) %*% W;
    trSigma_WX = sum( diag( Sigma_X %*% Phi ) );

    # restructure the linear term of the objective function 
    FirstDegree = ( -1 / trSigma_WX ) * matrix( t( Phi %*% Sigma_XF ), N*K, );

    # restructure the quadratic term of the objective function
    SecondDegree = Sigma_F;
    for( k in 1 : (N - 1) )
    {
        SecondDegree = blkdiag( SecondDegree, Sigma_F );
    }
    SecondDegree = ( 1 / trSigma_WX ) * t(kron( sqrt( Phi ), diag( 1, K ))) %*% SecondDegree %*% kron( sqrt(Phi), diag( 1, K ) );

    # restructure the equality constraints 
    if( !length(Aeq)  )
    {
        AEq = Aeq;
    }else
    {
        AEq = NULL;
        for( k in 1 : N )
        {
            AEq = cbind( AEq, kron( diag( 1, K ), Aeq[k] ) );
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
    G = matrix( attributes(b)$primal, N, ) ;

    return( G );
}