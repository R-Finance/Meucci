
#' Generates a multivariate i.i.d. sample of lenght J from the normal-inverse-Wishart distribution, as described in 
#' A. Meucci "Risk and Asset Allocation", Springer, 2005.
#'
#'  @param  Mu_0     : [vector]
#'  @param  T_0      : [scalar]
#'  @param  Sigma_0  : [matrix]
#'  @param  nu_0     : [scalar]
#'  @param  J        : [scalar]
#'  
#'  @return Mu       : [vector]
#'  @return Sigma    : [matrix]
#'  @return InvSigma : [matrix]
#'
#'  @note
#'  Mu|Sigma   ~ N(Mu_0,Sigma/T_0)
#'  inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "RandNormalInverseWishart.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

RandNormalInverseWishart = function(Mu_0, T_0, Sigma_0, nu_0, J)
{
  N = length( Mu_0 );
  VecIndex = NULL;
  for( n in 1 : N )
  {
      VecIndex[ n ] = cbind( VecIndex, ( n-1 ) * N +( n:N ) ); ##ok<AGROW>
  }

  invSigma_0 = solve(Sigma_0) %*% diag( 1, dim( Sigma_0 ));
  Phi = invSigma_0 / nu_0;

  Mu       = NULL;
  Sigma    = NULL;
  InvSigma = NULL;

  for( j in 1 : J )
  { 
      Inv_Sigma = rwishart( df = nu_0, Sigma = Phi );
      InvSigma  = rbind( InvSigma, Inv_Sigma[ VecIndex ] );    
               
      S = solve(Inv_Sigma) %*% diag( 1, dim( Inv_Sigma ) );
      Sigma = rbind( Sigma, S[VecIndex] );
      
      M = rmvnorm( nrow(Mu_0) * nrow(S/T_0), Mu_0, S/T_0);
      Mu = rbind( Mu, M );
  }

  return( list( Mu = Mu, Sigma = Sigma , InvSigma = InvSigma  ) )
}

