#' This script displays the cdf of the copula of a normal distribution, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}, 
#' "E 35 - Cdf of the normal copula".
#'
#' See Meucci's script for "S_DisplayNormalCopulaCdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' 

#############################################################################################################
### Input parameters
Mu = c( 0, 0 );     
r  = -0.999;            
sigmas = c(1, 1 );    
Sigma = diag( sigmas ) %*% rbind( c( 1, r ), c( r, 1 ) ) %*% diag( sigmas );

#############################################################################################################
### Grid
GridSide1 = seq( 0.05, 0.95, 0.05 );
GridSide2 = GridSide1;
nMesh = length(GridSide1);

#############################################################################################################
### Compute cdf of copula

F_U = matrix( NaN, nMesh, nMesh);

for ( j in 1 : nMesh )
{
    for ( k in 1 : nMesh)
    {
        u = c( GridSide1[ j ], GridSide2[ k ] );
        x= qnorm( u, Mu, sigmas );        
        F_U[ j, k ] = pmvnorm( lower = -Inf, upper = x, mean = Mu, corr = Sigma );         
    }
}

#mesh representation    

persp( GridSide1, GridSide2, F_U,
	theta = 7 * 45, phi = 30, expand=0.6, col='lightblue', shade=0.75, ltheta=120, 
	ticktype='detailed', xlab = "U_1", ylab = "U_2", zlab = "copula cdf" );