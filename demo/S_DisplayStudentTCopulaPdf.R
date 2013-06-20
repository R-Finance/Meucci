#'This script displays the pdf of the copula of a Student t distribution, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_DisplayStudentTCopulaPdf.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#############################################################################################################
### input parameters

Mu = rbind( 0, 0 );     
r  = 0.5;            
sigmas = rbind( 1, 2 );    
Sigma = diag( c( sigmas ) ) %*% rbind( c( 1, r ), c( r, 1 ) ) %*% diag( c( sigmas ) );
#nu = 1; Sigma(1,2) = 0; Sigma(2,1) = 0;
nu = 200; 

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
        f_U[ j, k ] = StudentTCopulaPdf( u, nu, Mu, Sigma );         
    }
}

#mesh representation    

persp( GridSide1, GridSide2, f_U,
	theta = 7 * 45, phi = 30, expand=0.6, col='lightblue', shade=0.75, ltheta=120, 
	ticktype='detailed', xlab = "U_1", ylab = "U_2", zlab = "copula pdf" );
