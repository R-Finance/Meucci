#' This script computes the multivariate shrinkage estimators of location and scatter under the normal assumption,
#'  as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_ShrinkageEstimators.m"
#'

##################################################################################################################
### Inputs
N  = 5;
T  = 30;
Mu = runif(N);
A  = matrix(runif( N* N), N, N) - 0.5;
Sigma = A %*% t(A);
    
#################################################################################################################
### Generate normal sample
X = rmvnorm( T, Mu, Sigma );

# estimate sample parameters
Mu_hat = matrix( apply( X, 2, mean ));
Sigma_hat = cov( X ) * (T - 1) / T;

#################################################################################################################
### Shrinkage of location parameter 

# target 
b = matrix( 0, N, 1 ); 

# compute optimal weight
Lambda_hat = eigen( Sigma_hat )$values; 
a = 1 / T * ( sum( Lambda_hat ) - 2 * max( Lambda_hat) ) / ( t( Mu_hat-b)  %*% ( Mu_hat - b ) ); 

# restrict to sensible weight
a = max( 0, min( a, 1) );     

# shrink
Mu_shr = ( 1 - a ) * Mu_hat + a * b;

print(Mu_hat);
print(Mu_shr);

#################################################################################################################
### Shrinkage of scatter parameter

# target
C = mean(Lambda_hat) * diag( 1, N );

# compute optimal weight
Numerator = 0;
for( t in 1 : T )
{
    Numerator = Numerator + 1 / T * sum(diag( (matrix(X[ t, ])%*%X[ t, ] - Sigma_hat ) ^ 2 ));
}

Denominator = sum( diag( ( Sigma_hat - C ) ^ 2 ));
a = 1 / T * Numerator / Denominator;

# restrict to sensible weight
a = max(0, min(a, 1)); 

# shrink
Sigma_shr = ( 1 - a ) * Sigma_hat + a * C;

print(Sigma_hat);
print(Sigma_shr);

