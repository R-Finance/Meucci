#' Generates normal simulations whose sample moments match the population moments
#'
#' Adapted from file 'MvnRnd.m'. Most recent version of article and code available at http://www.symmys.com/node/162
#' see A. Meucci - "Simulations with Exact Means and Covariances", Risk, July 2009
#'
#' @param M         a numeric indicating the sample first moment of the distribution
#' @param S         a covariance matrix
#' @param J         a numeric indicating the number of trials
#'
#' @author Ram Ahluwalia \email{rahluwalia@@gmail.com}
#' @references 
#' \url{http://www.symmys.com}
#' TODO: Add Schur decomposition. Right now function is only sampling from mvrnorm so sample moments do no match population moments
#' I have sample code commented out below to implement this correctly but I require a function that returns the unitaryMatrix from a Schur decomposition
#' @export
MvnRnd = function( M , S , J )
{
  library(MASS)
  X = MASS::mvrnorm( n = J , mu = M , Sigma = S ) # Todo: need to swap with Meucci function and Schur method
  return( X = X )
    
    # # compute sample covariance: NOTE defined as "cov(Y,1)", not as "cov(Y)"
    # S_ = cov( Y , 1 )
    # 
    # # solve Riccati equation using Schur method
    #     zerosMatrix = matrix( rep( 0 , length( N * N ) ) , nrow = N )
    #     # define the Hamiltonian matrix
    #     H1 = cbind( zerosMatrix , -1*S_ )
    #     H2 = cbind( -S , zerosMatrix ) 
    #     H = rbind( H1 , H2 )
    #     # perform its Schur decomposition. 
    #     # TODO: check that the result returns real eigenvalues on the diagonal. ?Schur seems to give an example with real eigenvalues
    #     schurDecomp = Schur( H )
    #     T = SchurDecomp
    #     # U_ = unitaryMatrix??? TODO: Find a function in R that returns the unitaryMatrix from a Schur decomposition
    #     # U = ordschur(U_,T_,'lhp')
    #     # U_lu = U(1:N,1:N)
    #     # U_ld = U(N+1:end,1:N)
    #     # B = U_ld/U_lu
    # 
    # # affine transformation to match mean and covariances
    # # X = Y%*%B + repmat(M', J , 1 ) 
    #
}