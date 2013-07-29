library( mvtnorm );

#' Simulate invariants for the regression model, as described in  A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005, Chapter 4.
#'
#'  @param   mu_x    : [scalar] 
#'  @param   sig_x   : [scalar] std for x
#'  @param   nu_f    : [scalar] dof for x
#'  @param   sig_f   : [scalar] std for x
#'  @param   nu_w    : [scalar] dof for w
#'  @param   Sigma_w : [matrix] ( 2 x 2 )  covariance matrix 
#'  @param   nu_w    : [scalar] dof for w
#'  
#'  @return  X       : [vector] ( J x 1 ) 
#'  @return  F       : [vector] ( J x 1 ) 
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "GenerateInvariants.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

GenerateInvariants = function( mu_x, sig_x, nu_f, sig_f, nu_w, Sigma_w, J ) 
{
    diag_W = 0;

    for( s in 1 : nu_w ) 
    {
        Z = rmvnorm( J, rbind( 0, 0 ), Sigma_w );
        diag_W = diag_W + Z * Z;
    }

    a_w = nu_w / 2;
    b_w = 2 * diag( Sigma_w );
    a_f = nu_f / 2;
    b_f = 2 * sig_f ^ 2;
    U_x = pgamma( diag_W[ , 1 ], a_w, b_w[ 1 ] );
    X   = qlnorm( U_x, mu_x, sig_x );
    U_f = pgamma( diag_W[ , 2 ], shape = a_w, scale = b_w[ 2 ] );
    F   = qgamma( U_f, shape = a_f, scale = b_f );

    return( list( X = matrix(X), F = matrix(F) ) );
}


#' This script simulates statistics for a regression model and compare it theoretical ones, 
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_TStatApprox.m"
#'

##################################################################################################################
### Inputs
T = 25;
J = 1500;
mu_x  = 0.1 * ( runif(1)  - 0.5 );
sig_x = 0.2 * runif(1);
nu_f  = ceiling( 10 * runif(1) );
sig_f = 0.2 * runif(1);
nu_w  = ceiling( 10 * runif(1) );
dd    = matrix( runif(4), 2, 2 )  - 0.5;
Sigma_w = dd %*% t(dd);

##################################################################################################################
### Compute market features in simulation
GI    = GenerateInvariants( mu_x, sig_x, nu_f, sig_f, nu_w, Sigma_w, J );

Mu    = apply( cbind(GI$X, GI$F), 2, mean );
Sigma = cov( cbind( GI$X, GI$F ) );
mu_X  = Mu[ 1 ];
mu_F  = Mu[ 2 ];
sig_X = sqrt( Sigma[ 1, 1 ] );
sig_F = sqrt( Sigma[ 2, 2 ] );
rho = Sigma[ 1, 2 ]  / sqrt( Sigma[ 1, 1 ]  * Sigma[ 2, 2 ] );

Alpha = mu_X - mu_F * rho * sig_X / sig_F;
Beta  = rho * sig_X / sig_F;
sig   = sig_X * sqrt( 1 - rho^2 );

##################################################################################################################
### Randomize time series and compute statistics as random variables
mu_X_hat   = matrix( 0, 1, J);
sig2_X_hat = matrix( 0, 1, J);
Alpha_hat  = matrix( 0, 1, J);
Beta_hat   = matrix( 0, 1, J);
sig2_a_hat = matrix( 0, 1, J);
sig2_b_hat = matrix( 0, 1, J);
sig2_U_hat = matrix( 0, 1, J);
Sigma_F    = matrix( 0, 2, 2);
for( j in 1 : J )
{
    GI = GenerateInvariants( mu_x, sig_x, nu_f, sig_f, nu_w, Sigma_w, T );
    
    # t-stat for mean
    mu_X_hat[ j ]    = mean( GI$X );
    sig2_X_hat[ j ]  = (dim(GI$X)[1]-1)/dim(GI$X)[1]* var( GI$X);

    # t-stat for regression
    Sigma_XF = cbind( mean( GI$X ),  mean( GI$X * GI$F ) );
    Sigma_F[ 1, 1 ]  = 1;
    Sigma_F[ 1, 2 ]  = mean( GI$F );
    Sigma_F[ 2, 1 ]  = Sigma_F[ 1, 2 ];
    Sigma_F[ 2, 2 ]  = mean( GI$F * GI$F );
    inv_Sigma_F      = solve(Sigma_F);
    sig2_a_hat[ j ]  = inv_Sigma_F[ 1, 1 ];
    sig2_b_hat[ j ]  = inv_Sigma_F[ 2, 2 ];

    AB_hat           = Sigma_XF %*% inv_Sigma_F;
    Alpha_hat[ j ]   = AB_hat[ 1 ];
    Beta_hat[ j ]    = AB_hat[ 2 ];        
    
    X_ = Alpha_hat[ j ] + Beta_hat[ j ] * GI$F;
    U  = GI$X - X_;
    sig2_U_hat[ j ]  = (dim(U)[1]-1)/dim(U)[1] * var( U ); #MOMENT
}

t_m = ( mu_X_hat - mu_X ) / sqrt( sig2_X_hat / ( T - 1 ) );
t_a = ( Alpha_hat - Alpha ) / sqrt( sig2_a_hat * sig2_U_hat / ( T - 2 ) );
t_b = ( Beta_hat - Beta ) / sqrt( sig2_b_hat * sig2_U_hat / ( T - 2 ) );

##################################################################################################################
### Display results
# should be uniform
dev.new();
par( mfrow = c( 3, 1) );

NumBins = round( 10 * log( J ) );

U_m     = pt( t_m, T-1 );
hist( U_m, NumBins );

U_a= pt( t_a, T-2 );
hist( U_a, NumBins );

U_b = pt( t_b, T-2 );
hist( U_b, NumBins );

# qqplots for comparison
dev.new();
par( mfrow = c( 1, 3 ) )
qqplot( t_m, matrix( rt( length( t_m ), T-1 ), dim( t_m ) ) );
qqplot( t_a, matrix( rt( T-2, length( t_a ) ), dim( t_a ) ) );
qqplot( t_b, matrix( rt( T-2, length( t_b ) ), dim( t_b ) ) );

