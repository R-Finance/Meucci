#' Uses Entropy Pooling to compute a double-decay covariance matrix.
#'
#' This function uses Entropy Pooling to compute a double-decay covariance matrix, as described in  
#' A. Meucci, "Personalized Risk Management: Historical Scenarios with Fully Flexible Probabilities"
#' GARP Risk Professional, Dec 2010, p 47-51
#' 
#' @param   X   matrix representing the risk drivers.
#' @param   m   matrix of zeros, representing the expectation of the risk drivers.
#' @param   S   matrix representing the double-decay estimation for the correlation matrix of the risk drivers.
#' @return  p   list containing the vector of posterior probabilities and information about the optimization performance.
#' 
#' @references 
#' \url{http://www.symmys.com/node/150}
#' See Meucci script for "S_MainFullFlexProbs.m"
#' 
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

Fit2Moms = function( X, m, S)
{
	N = dim(X);

	Aeq = matrix( 1, 1, N[1] );  # constrain probabilities to sum to one...
	beq = 1;

	Aeq = rbind( Aeq , t(X) ); # ...constrain the first moments...
	beq = rbind( beq, m );

	SecMom = S + m %*% t(m);  #...constrain the second moments...

	for ( k in  1:N[2] )
	{
		for ( l in k:N[2] )
		{
			Aeq = rbind( Aeq , t(X[ ,k] * X[ ,l] ) );
			beq = rbind( beq, SecMom[k,l] );
		}
	}

	p_0 = matrix( 1, N[1], 1) / N[1];

	return ( p = EntropyProg( p_0, matrix( , 0, 0), matrix( , 0, 0), Aeq , beq)$p_); # ...compute posterior probabilities

}