#' This script fits a time-series linear factor computing the industry factors loadings, as described in A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 114 – Time series factors: unconstrained time series industry factors".
#'
#' See Meucci's script for "S_TimeSeriesIndustries.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Loads weekly stock returns X and indices stock returns F
data("securitiesTS");
Data_Securities = securitiesTS$data[ , -1 ]; # 1st column is date

data("sectorsTS");
Data_Sectors = sectorsTS$data[ , -(1:2) ]; #1st column is for date, 2nd column is SPX index

##################################################################################################################
### Estimation
# linear returns for stocks
X = diff( Data_Securities ) / Data_Securities[ -nrow(Data_Securities), ];

# linear return for the factors
F = diff(Data_Sectors) / Data_Sectors[ -nrow(Data_Sectors), ];

T = dim(X)[1];
N = dim(X)[2];
K = dim(F)[ 2 ];

# compute sample estimates
E_X = matrix( apply(X, 2, mean ) );
E_F = matrix( apply(F, 2, mean ) );

XF = cbind( X, F )
SigmaJoint_XF = (dim(XF)[1]-1)/dim(XF)[1] * cov(XF);
Sigma_X       = SigmaJoint_XF[ 1:N, 1:N ];
Sigma_XF      = SigmaJoint_XF[ 1:N, (N+1):(N+K) ];
Sigma_F       = SigmaJoint_XF[ (N+1):(N+K), (N+1):(N+K) ];
Corr_F  	  = cov2cor(Sigma_F);
Corr_F  	  = tril(Corr_F, -1);

# compute OLS loadings for the linear return model
X_ = X - repmat( t(E_X), T, 1 );
F_ = F - repmat( t(E_F), T, 1 );
B  = Sigma_XF %*% solve(Sigma_F);
U  = X_ - F_ %*% t(B);

##################################################################################################################
### Residual analysis

UF = cbind(U,F);
SigmaJoint_UF = ( dim( UF )[1]-1 )/dim( UF )[1] * cov( UF );
CorrJoint_UF  = cov2cor( SigmaJoint_UF );

# correlations of residuals with factors is null
Corr_UF      = CorrJoint_UF[ 1:N, (N+1):(N+K) ];
mean_Corr_UF = mean( abs( as.array( Corr_UF )) );
max_Corr_UF  = max( abs( as.array( Corr_UF )) );

disp( mean_Corr_UF );
disp( max_Corr_UF );

dev.new();
hist( Corr_UF, 100);

# correlations between residuals is not null
Corr_U = tril( CorrJoint_UF[ 1:N, 1:N ], -1 );
Corr_U = Corr_U[ Corr_U != 0 ];
mean_Corr_U = mean( abs( Corr_U ) );
max_Corr_U  = max( abs( Corr_U ) );
disp( mean_Corr_U );
disp( max_Corr_U );

dev.new();
hist(Corr_U, 100);

