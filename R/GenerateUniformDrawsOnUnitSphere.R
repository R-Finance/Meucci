#' Generate a uniform sample on the unit hypersphere, as described in  A. Meucci,
#'  "Risk and Asset Allocation", Springer, 2005.
#'  
#'	@param   J : [scalar] number of draws
#'	@param   N : [scalar] dimension
#'  
#'	@return   X  : [matrix] (T x N) of draws
#'
#'@note
#' Initial script by Xiaoyu Wang - Dec 2006
#' We decompose X=U*R, where U is a uniform distribution on unit sphere and
#     R is a distribution on (0,1) proportional to r^(Dims-1), i.e. the area of surface of radius r 
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "GenerateUniformDrawsOnUnitSphere.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

GenerateUniformDrawsOnUnitSphere = function(J, N)
{
l = matrix( 1, 1, N );

# 1. Generate U
Z = matrix( runif(J*N), J, N );
normZ = sqrt( apply( Z * Z, 1, sum ) );
U = Z / ( normZ %*% l );

# 2. Generate R
# the pdf of R is proportional to r^(N-1) therefore the cdf of R is r^N
# we use quantile function of R sample R from uniform simulations
Y = runif(J);
R = Y ^ ( 1/N );

# 3. Generate X
X = U * ( R %*% l );
	return( X );
}