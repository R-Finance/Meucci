#' This script fits a cross-sectional linear factor model creating industry factors, as described in A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_CrossSectionIndustries.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data
# loads weekly stock returns X and indices stock returns F
load(" ../data/securitiesTS.Rda");
Data_Securities = securitiesTS$data[ , -1 ]; # 1st column is date

load("../data/securitiesIndustryClassification.Rda");
Securities_IndustryClassification = securitiesIndustryClassification$data;

##################################################################################################################
### Estimation
# linear returns for stocks
X = diff( Data_Securities ) / Data_Securities[ -nrow(Data_Securities), ];

T = dim(X)[1];
N = dim(X)[2];
B = Securities_IndustryClassification[ 1:N, ];
K = dim(B)[ 2 ];

# compute sample estimates
E_X = matrix( apply(X, 2, mean ) );
Sigma_X = (dim(X)[1]-1)/dim(X)[1] * cov(X);

# The optimal loadings turn out to be the standard multivariate weighted-OLS.
Phi = diag(1 / diag( Sigma_X ), length(diag( Sigma_X ) ) );
tmp = t(B) %*% Phi %*% B ;
F   = t( diag( diag(tmp) ^( -1 ), dim(tmp) ) %*% t(B) %*% Phi %*% t(X));

# compute intercept a and residual U
E_F = matrix( apply( F, 2, mean ) );
a   = E_X - B %*% E_F;
A_  = repmat( t(a), T, 1);
U   = X - A_ - F %*% t(B);

##################################################################################################################
### Residual analysis
M = cbind( U, F );
SigmaJoint_UF = (dim(M)[1]-1)/dim(M)[1] * cov( M );
CorrJoint_UF  = cov2cor(SigmaJoint_UF);
Sigma_F = (dim(F)[1]-1)/dim(F)[1] * cov(F);
Corr_F  = cov2cor( Sigma_F );
Corr_F  = tril(Corr_F, -1);
Corr_F  = Corr_F[ Corr_F != 0 ];

Corr_U = tril(CorrJoint_UF[ 1:N, 1:N ], -1);
Corr_U = Corr_U[ Corr_U != 0 ];
mean_Corr_U = mean( abs(Corr_U));
max_Corr_U  = max( abs(Corr_U));
disp(mean_Corr_U);
disp(max_Corr_U);

dev.new();
hist(Corr_U, 100);

Corr_UF = CorrJoint_UF[ 1:N, (N+1):(N+K) ];
mean_Corr_UF = mean( abs(Corr_UF ) );
max_Corr_UF  = max( abs(Corr_UF ) );
disp(mean_Corr_U);
disp(max_Corr_U);

dev.new();
hist(Corr_UF, 100);

