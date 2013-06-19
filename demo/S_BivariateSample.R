#' This script generates draws from a bivariate distribution with different marginals,
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' \url{http://}
#' See Meucci's script for "S_AnalyzeLognormalCorrelation.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

library(mvtnorm);
library(latticeExtra);

###################################################################################################################
### input parameters

nSim = 10000;

# input for bivariate normal distribution

NormCorr   = -0.8;
NormStDev  = rbind( 1, 3 );  # NOTE: this input plays no role in the final output
NormExpVal = rbind( -2, 5 ); # NOTE: this input plays no role in the final output

# input for first marginal
nu_1 = 9;
sigmasq_1 = 2;

mu_2 = 0;
sigmasq_2 = 0.04;

# input for second marginal
nu_2 = 7;

###################################################################################################################
### Generate draws from a bivariate normal distribution

NormCorrMatrix = rbind( c( 1, NormCorr ), c( NormCorr, 1 ));
NormCovMatrix  = diag( c( NormStDev ) ) %*% NormCorrMatrix %*% diag( c( NormStDev) );

Z = rvnorm( nSim, NormExpVal, NormCovMatrix );

Z_1 = Z[, 1];
Z_2 = Z[, 2];

# display marginals: as expected, they are normal

NumBins = round(10 * log(nSim));
par( mfrow = c( 2, 1) );
hist( Z_1, NumBins, xlab = "normal 1", ylab = "" );
hist( Z_2, NumBins, xlab = "normal 2", ylab = "" );

plot( Z_1, Z_2, type = "p", xlab = "normal 1", ylab = "normal 2" );


# 3d histograms 

NumBins2D = round(sqrt(100 * log(nSim)));
Z_3 = table( cut (Z_1, NumBins2D ), cut ( Z_2, cloud ));

cloud( Z_3, panel.3d.cloud = panel.3dbars, scales = list( arrows = FALSE, just = "right" ), 
	xlab = "normal 1", ylab = "normal 2", zlab="", main = "pdf normal" );

###################################################################################################################
### Generate draws from the copula

U_1 = pnorm( Z[ , 1 ], NormExpVal[ 1 ], NormStDev[ 1 ]);  # grade 1
U_2 = pnorm( Z[ , 2 ], NormExpVal[ 2 ], NormStDev[ 2 ]);  # grade 2
U = c( U_1, U_2 ); # joint realizations from the required copula

# plot copula
NumBins = round(10 * log(nSim));
par( mfrow = c( 2, 1) );
hist( U_1, NumBins, xlab = "grade 1", ylab = "", main = "" );
hist( U_2, NumBins, xlab = "grade 2", ylab = "", main = "" );

# joint sample
plot(U_1, U_2, xlab="grade 1", ylab="grade 2" );

# 3d histogram
NumBins2D = round(sqrt(100 * log(nSim)));
U_3 = table( cut (U_1, NumBins2D ), cut ( U_2, NumBins2D ));
cloud( U_3, panel.3d.cloud = panel.3dbars, scales = list( arrows = FALSE, just = "right" ), 
	xlab = "grade 1", ylab = "grade 2", zlab="", main = "pdf copula" );

###################################################################################################################
### Generate draws from the joint distribution
a = nu_1 / 2;
b = 2 * sigmasq_1;
X_1 = qgamma( U_1, a, b );

sigma_2 = sqrt( sigmasq_2 );
X_2 = qlnorm( U_2, mu_2, sigma_2 );

X = C(X_1, X_2); # joint realizations from the required distribution

###################################################################################################################
### Plot joint distribution
# marginals: as expected, the histograms (pdf's) do NOT change as NormCorr varies

NumBins = round(10 * log(nSim));


par( mfrow = c( 2, 1) );
# Student t distribution
hist( X_1, NumBins, xlab = "gamma", ylab = "", main = "" );
# chi-square distribution
hist( X_2, NumBins, xlab = "lognormal", ylab = "", main = "" );

# joint sample
plot(X_1, X_2, xlab="gamma", ylab="lognormal" );

# 3d histogram
NumBins2D = round(sqrt(100 * log(nSim)));
X_3 = table( cut (X_1, NumBins2D ), cut ( X_2, NumBins2D ));
cloud( X_3, panel.3d.cloud = panel.3dbars, scales = list( arrows = FALSE, just = "right" ), 
	xlab = "gamma", ylab = "lognormal", zlab="", main = "pdf joint distribution" );