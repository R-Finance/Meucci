
#'This script displays the pdf of the copula of a lognormal distribution, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_DisplayLognormalCopulaPdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#############################################################################################################
### Input parameters
Mu = rbind( 100.0, -30.0 );     
r  = 0.8;            
sigmas = rbind( 1000, 0.01 );    
nu = 100;
Sigma = diag( c( sigmas ) ) %*% rbind( c( 1, r ), c( r, 1 ) ) %*% diag( c( sigmas ) );

#############################################################################################################
### Grid
GridSide1 = seq( 0.05, 0.95, 0.05 );
GridSide2 = GridSide1;
nMesh = length(GridSide1);

#############################################################################################################
### Compute pdf of copula

f_U = matrix( NaN, nMesh, nMesh);

for ( j in 1 : nMesh )
{
    for ( k in 1 : nMesh)
    {
        u = c( GridSide1[ j ], GridSide2[ k ] );        
        f_U[ j, k ] = LognormalCopulaPdf(u, Mu, Sigma);         
    }
}

#mesh representation    

persp( GridSide1, GridSide2, f_U,
	theta = 7 * 45, phi = 30, expand=0.6, col='lightblue', shade=0.75, ltheta=120, 
	ticktype='detailed', xlab = "U_1", ylab = "U_2", zlab = "copula pdf" );
