#' This script familiarizes the users with multivariate Bayesian estimation.
#' A normal time series is generated a normal-inverse-Wishart prior is set.
#' The ensuing normal-inverse-Wishart posterior is computed and analyzed numerically and analytically.
#' Described in A. Meucci,"Risk and Asset Allocation",Springer, 2005, Chapter 7.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 282 - Bayesian: normal-inverse-Wishart posterior".
#'
#' See Meucci's script for "S_AnalyzeNormalInverseWishart.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

###################################################################################################################
### Inputs
N = 1; # market dimension
nSim = 2000; 

###################################################################################################################
### History : X_t ~ N(M,S), t=1,...,T
# set parameters
M = 1 * array( 1, N );
s = 1 * array( 1, N );
r = 0.7;
C = (1 - r) * diag( 1, N ) + r * matrix( 1, N, N);
S = diag( s, N ) %*% C %*% diag( s, N );
T = 520;

# generate time series
X = rmvnorm( T, M, S);
Mu_Hat = apply(X, 2, mean);
Sigma_Hat = cov(X);

###################################################################################################################
### Prior : Mu|Sigma ~ N(Mu_0,Sigma/T_0)
###         Omega == inv(Sigma) ~ W(Nu_0,inv(Sigma_0)/Nu_0)
# set parameters
Mu_0 = 2 * array( 1, N );
T_0  = 520;
s_0  = 2 * array( 1, N );
r    = 0.3;
C    = ( 1 - r ) * diag( 1, N ) + r * matrix( 1, N, N );
Sigma_0 = diag( s_0, N ) %*% C %*% diag( s_0, N );
Nu_0 = 520;

# generate simulations
RNIWPrior= RandNormalInverseWishart(Mu_0, T_0, Sigma_0, Nu_0, nSim);

# plot results
PlotMarginalsNormalInverseWishart( RNIWPrior$Mu , RNIWPrior$InvSigma, Mu_0, T_0, Sigma_0, Nu_0, "prior" );

###################################################################################################################
### Posterior : Mu|Sigma ~ N(Mu_1,Sigma/T_1)
###             Omega == inv(Sigma) ~ W(Nu_1,inv(Sigma_1)/Nu_1)

# set parameters
T_1  = T_0 + T;
Mu_1 = ( T_0 %*% Mu_0 + T %*% Mu_Hat) / T_1;
Nu_1 = Nu_0 + T;
Sigma_1 = ( Nu_0 %*% Sigma_0 + T %*% Sigma_Hat + ( T %*% T_0 / (T + T_0)) %*% (Mu_0 - Mu_Hat) %*% t(Mu_0 - Mu_Hat) ) / Nu_1;

# generate simulations
RNIWPost = RandNormalInverseWishart(Mu_1, T_1, Sigma_1, Nu_1, nSim);

# plot results
PlotMarginalsNormalInverseWishart( RNIWPost$Mu, RNIWPost$InvSigma, Mu_1, T_1, Sigma_1, Nu_1, "posterior" );

###################################################################################################################
### Compute statistics
Mu_CE_Num = apply( RNIWPost$Mu, 2, mean);
Mu_CE_Anal = t( Mu_1 );
Mu_Hat = t( Mu_Hat );
Mu_0 = t( Mu_0 );

Mu_Scatter_Num  = cov( RNIWPost$Mu );
Mu_Scatter_Anal = Nu_1 / ( Nu_1 - 2 ) * Sigma_1 / T_1;

Sigma_CE_Num  = apply(RNIWPost$Sigma,2, mean);
Sigma_CE_Anal = Sigma_1;
print(Sigma_Hat);
print(Sigma_0);

Sigma_Scatter_Num = cov(RNIWPost$Sigma);

InvSigma_CE_Num = apply(RNIWPost$InvSigma, 2, mean);
S = solve( Sigma_1 );
InvSigma_CE_Anal = S;

