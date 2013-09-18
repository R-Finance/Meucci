#' This script fits a time-series linear factor computing the industry factors loadings,  where the loadings are 
#' bounded and constrained to yield unit exposure, as described in A. Meucci, "Risk and Asset Allocation",
#' Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 115 – Time series factors: generalized time-series industry factors".
#'
#' See Meucci's script for "S_TimeSeriesConstrainedIndustries.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Loads weekly stock returns X and indices stock returns F
data("securitiesTS");
Data_Securities = securitiesTS$data[ , -1 ]; # 1st column is date

data("sectorsTS");
Data_Sectors = sectorsTS$data[ , -(1:2) ]; #1st column is date, 2nd column is SPX

##################################################################################################################
### Estimation
# linear returns for stocks
X = diff( Data_Securities ) / Data_Securities[ -nrow(Data_Securities), ];
X = X[ ,1:20 ]; # consider first stocks only to lower run time (comment this line otherwise)

# linear return for the factors
F = diff(Data_Sectors) / Data_Sectors[ -nrow(Data_Sectors) , ];

T = dim(X)[1];
N = dim(X)[2];
K = dim(F)[ 2 ];

##################################################################################################################
# Solve for B that maximises sample r-square of F'*B' with X
# under constraints A*B<=D and Aeq*B=Deq (A,D, Aeq,Deq conformable matrices)
W   = diag( 1, N );
A 	= NULL;
D   = NULL;
Aeq = matrix( 1, K, N ) / N
Deq = matrix( 1, K, 1 );
lb  = 0.8
ub  = 1.2
B   = MaxRsqTS(X, F, W, A, D, Aeq, Deq, lb, ub);

# compute sample estimates
E_X = matrix( apply(X, 2, mean ) );
E_F = matrix( apply( F, 2, mean ) );
XF = cbind( X, F );
SigmaJoint_XF = (dim(XF)[1]-1)/dim(XF)[1] * cov(XF);
	
Sigma_X  = SigmaJoint_XF[ 1:N, 1:N ];
Sigma_XF = SigmaJoint_XF[ 1:N, (N+1):(N+K) ];
Sigma_F  = SigmaJoint_XF[ (N+1):(N+K), (N+1):(N+K) ];

# compute intercept a and residual U
a = E_X - B %*% E_F;
U = X - repmat(t(a), T, 1) - F %*% t(B);


##################################################################################################################
### Residual analysis

M = cbind( U, F );
SigmaJoint_UF = (dim(M)[1]-1)/dim(M)[1] * cov( M );
CorrJoint_UF  = cov2cor(SigmaJoint_UF);

# correlation between residuals is not null
Corr_U = tril(CorrJoint_UF[ 1:N, 1:N ], -1);
Corr_U = Corr_U[ Corr_U != 0 ];
mean_Corr_U = mean( abs(Corr_U));
max_Corr_U  = max( abs(Corr_U));

dev.new();
hist(Corr_U, 100);

# correlation between residuals and factors is not null
Corr_UF = CorrJoint_UF[ 1:N, (N+1):(N+K) ];
mean_Corr_UF = mean( abs(Corr_UF ) );
max_Corr_UF  = max( abs(Corr_UF ) );

dev.new();
hist(Corr_UF, 100);