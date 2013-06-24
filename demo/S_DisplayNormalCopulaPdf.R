#'This script displays the pdf of the copula of a normal distribution, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_DisplayNormalCopulaPdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#############################################################################################################
### input parameters
Mu = rbind( 1,  -1 );     
r  = 0.7;            
sigmas = c( 1, 1 );    
Sigma = diag( sigmas ) %*% rbind( c( 1, r ), c( r, 1 ) ) %*% diag( sigmas );

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
        f_U[ j, k ] = NormalCopulaPdf(u, Mu, Sigma);         
    }
}

#mesh representation    

persp( GridSide1, GridSide2, f_U,
	theta = 7 * 45, phi = 30, expand=0.6, col='lightblue', shade=0.75, ltheta=120, 
	ticktype='detailed', xlab = "U_1", ylab = "U_2", zlab = "copula pdf" );
