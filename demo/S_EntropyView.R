#' This script  illustrates the Entropy Pooling approach, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 9.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 304 - Entropy pooling".
#'
#' See Meucci's script for "S_EntropyView.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Market simulations
nSim = 100000;
B = ( runif( nSim ) < 0.5);
X = B * rnorm( nSim, -1, 1 ) + ( 1 - B ) * rnorm( nSim, 1, 1 );

##################################################################################################################
### View
# specify view E{X} = 0.5 and constraint 1'*p = 1.
p_prior = matrix( 1, nSim, 1) / nSim;
Aeq = rbind( X, matrix( 1, 1, nSim ) );
beq = rbind( 0.5, 1 );

##################################################################################################################
### Posterior market distribution using the Entropy Pooling approach
#Using package's EntropyProg instead of Books EntropyMinimization (Same function, different names)
p_post = EntropyProg( p_prior, Aeq = Aeq, beq = beq)$p;
pHistPriorPosterior(X,p_prior, p_post);
fprintf('prior sample mean = #f\n', mean(X));
fprintf('posterior sample mean = #f\n', X' * p_post);

