#' This script generates a sample from the 2x2 Wishart distribution.
#' it shows that determinant and trace are positive, i.e. the matrix is positive
#' it shows that the marginal diagonal are gamma-distributed
#' Described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_Wishart.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export


library(scatterplot3d);

###################################################################################################################
### Set inputs
s  = c( 1, 1 ); # variances
r  = 0.3; # correlation
Sigma = diag( c( s ) ) %*% rbind( c( 1, r ), c( r, 1 ) ) %*% diag( c( s ) );
nu = 5; # degrees of freedom
nSim = 10000;

###################################################################################################################
### Generate draws

# initialize storage vectors/matrices
W_xx  = matrix( NaN, nSim, 1 ); 
W_yy  = matrix( NaN, nSim, 1 ); 
W_xy  = matrix( NaN, nSim, 1 ); 
Vec_W = matrix( NaN, nSim, 4 ); 
Dets  = matrix( NaN, nSim, 1 ); 
Traces = matrix( NaN, nSim, 1 ); 

# generate draws and store elements of W, trace and determinant

for (j in 1 : nSim )
{
  X = rmvnorm( nu, matrix( 0, 2, 1 ), Sigma );
  W = t(X) %*% X;
  Dets[ j ] = det( W );
  Traces[ j ] = sum(diag( W ));

  W_xx[ j ] = W[ 1, 1 ];
  W_yy[ j ] = W[ 2, 2 ];
  W_xy[ j ] = W[ 1, 2 ];

  Vec_W [ j, ] = as.vector( W );
}



#####################################################################################################################
### Check positivity of trace and determinant

NumBins = round(10 * log(nSim));
dev.new();
par( mfrow = c( 2, 1) );
hist(Traces, NumBins, xlab = "trace", ylab = "", main = "" );
hist(Dets, NumBins, xlab = "determinant", ylab = "", main = "" );

########################################################################################################################################
### Plot cloud of draws
# tri-variate joint

dev.new();
scatterplot3d(W_xx, W_yy, W_xy, xlab = expression( "W"[ 11 ] ),
 ylab = expression( "W"[ 22 ] ), zlab = expression( "W"[ 12 ] ) );

# bi-variate marginals

dev.new(); 
par( mfrow = c( 2, 2) );

plot(W_xx, W_xy, xlab = expression( "W"[ 11 ] ), ylab = expression( "W"[ 12 ] ) );

plot(W_yy, W_xy, xlab = expression( "W"[ 22 ] ), ylab = expression( "W"[ 12 ] ) );

plot(W_xx, W_yy, xlab = expression( "W"[ 11 ] ), ylab = expression( "W"[ 22 ] ) );

########################################################################################################################################
### Plot individual marginals histograms and pdfs

dev.new();
par( mfrow = c( 2, 2) );

D = hist(W_xx, NumBins, xlab = expression( "W"[ 11 ] ), ylab = "", main = "" );
RescalePdf = (0.5) * nSim;
y = dgamma( W_xx, shape = nu / 2, scale = 2 * Sigma[ 1, 1 ]) %*% RescalePdf;
points( W_xx, y, col = "blue"); 


hist(W_yy, NumBins, xlab = expression( "W"[ 22 ] ), ylab = "", main = "" )
RescalePdf = (D(2) - D(1)) * nSim;
y = dgamma(W_yy, shape = nu / 2, scale = 2 * Sigma[ 2, 2 ]) * RescalePdf;
points( W_yy, y, col = "blue"); 



hist(W_xy, NumBins, xlab = expression( "W"[ 12 ] ), ylab = "", main = "" );

#####################################################################################################################
### Compute summary statistics (analytical and sample-based)

Expected_Value = as.array(nu * Sigma);
K_22 = rbind(c( 1, 0,	0, 0 ), c( 0, 0, 1, 0 ), c( 0, 1, 0, 0 ), c( 0, 0, 0, 1 ) );

Covariance = nu * ( diag(1, 4) + K_22 ) %*% kronecker(Sigma, Sigma);

Sample_Mean = mean(Vec_W);
Sample_Covariance= cov(Vec_W);

print(Expected_Value);
print(Covariance);
print(Sample_Mean);
print(Sample_Covariance);

