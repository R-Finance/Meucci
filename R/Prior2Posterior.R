
#' Calculate the full-confidence posterior distributions of Mu and Sigma
#'
#' \deqn{ \tilde{ \mu }  \equiv \mu +  \Sigma  Q'    {\big(Q \Sigma  Q' \big)}^{-1}   \big( \tilde{\mu}_{Q} - Q \mu \big),
#' \\ \tilde{ \Sigma } \equiv \Sigma + \Sigma G' \big({\big(G \Sigma  G' \big)}^{-1} \tilde{ \Sigma }_G {\big(G \Sigma  G' \big)}^{-1} - {\big(G \Sigma  G' \big)}^{-1} \big) G \Sigma }
#' @param M     a numeric vector with the Mu of the normal reference model
#' @param Q     a numeric vector used to construct a view on expectation of the linear combination QX
#' @param M_Q   a numeric vector with the view of the expectations of QX
#' @param S     a covariance matrix for the normal reference model
#' @param G     a numeric vector used to construct a view on covariance of the linear combination GX
#' @param S_G   a numeric with the expectation associated with the covariance of the linear combination GX
#'
#' @return a list with 
#' @return M_   a numeric vector with the full-confidence posterior distribution of Mu
#' @return S_   a covariance matrix with the full-confidence posterior distribution of Sigma
#'
#' @references 
#' \url{http://www.symmys.com}
#' \url{http://ssrn.com/abstract=1213325}
#' A. Meucci - "Fully Flexible Views: Theory and Practice". See formula (21) and (22) on page 7
#' See Meucci script Prior2Posterior.m attached to Entropy Pooling Paper
#' @author Ram Ahluwalia \email{ram@@wingedfootcapital.com}
#' @export

Prior2Posterior = function( M , Q , M_Q , S , G , S_G )
{
  # Compute posterior moments
  
  if ( Q != 0 ) { M_ = M + S %*% t(Q) %*% solve( Q %*% S %*% t(Q) ) %*% ( M_Q - Q %*% M) }
  else { M_ = M }
  
  if ( G != 0 ) { S_ = S + (S %*% t(G)) %*% ( solve(G %*% S %*% t(G)) %*% S_G %*% solve(G %*% S %*% t(G)) - solve( G %*% S %*% t(G)) ) %*% (G %*% S) }
  else { S_ = S }
  
  return( list( M_ = M_ , S_ = S_ ) )
}