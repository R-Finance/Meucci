#' This script fits a cross-sectional linear factor model creating industry factors, where the industry factors 
#' are constrained to be uncorrelated with the market as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 3.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_CrossSectionConstrainedIndustries.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Loads weekly stock returns X and indices stock returns F
load("../data/securitiesTS.rda");
Data_Securities = securitiesTS$data[ , -1 ]; # 1st column is date

load("../data/securitiesIndustryClassification.rda");
Securities_IndustryClassification = securitiesIndustryClassification$data;
##################################################################################################################
### Linear returns for stocks
X = diff( Data_Securities ) / Data_Securities[ -nrow(Data_Securities), ];
X = X[ ,1:30 ]; # consider first stocks only to lower run time (comment this line otherwise)

T = dim(X)[1];
N = dim(X)[2];
B = Securities_IndustryClassification[ 1:N, ];
K = dim(B)[ 2 ];
m = array( 1/N, N );

# compute sample estimates
E_X = matrix( apply(X, 2, mean ) );
Sigma_X = (dim(X)[1]-1)/dim(X)[1] * cov(X);

# The optimal loadings turn out to be the standard multivariate weighted-OLS.
Phi = diag(1 / diag( Sigma_X ), length(diag( Sigma_X ) ) );
W 	= sqrt( Phi );

##################################################################################################################
### Solve for G that maximises sample r-square of X*G'*B' with X
###  under constraints A*G<=D and Aeq*G=Deq (A,D, Aeq,Deq conformable matrices)
A 	= NULL;
Aeq = t(m) %*% t(Sigma_X);
Deq = matrix( 0, K, 1 );
#BOUNDARIES 
lb = -100
ub = 700

#THE ROWS 3 and 6 should be 0 and instead of it we got outliers.
G   = MaxRsqCS(X, B, W, A, D, Aeq, Deq, lb, ub);

# compute intercept a and residual U
F   = X %*% t(G);
E_F = matrix(apply(F, 2, mean));
a   = E_X - B %*% E_F;
A_  = repmat(t(a), T, 1);
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