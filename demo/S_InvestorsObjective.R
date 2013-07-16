library(mvtnorm);
#' This script familiarizes the users with the objectives of different investors in a highly 
#' non-normal bi-variate  market of securities, as described in A. Meucci,"Risk and Asset 
#' Allocation",Springer, 2005,  Chapter 5.  
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_InvestorsObjective.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Parameters of first marginal
nu_1 = 3;
s2_1 = 3;

# parameters of second marginal
mu_2 = 0.1;
s2_2 = 0.2;

# correlation in normal copula
r = 0.5;

# number of simulations
J = 10000;  

# portfolio allocation
a = matrix(c( 1, 2 ));

# benchmark allocation
b = matrix( c( 2, 1 )); 

##################################################################################################################
### Compute current prices
p_1 = nu_1 * s2_1;
p_2 = exp( mu_2 + 0.5 * s2_2 ^ 2 );
p = matrix( c( p_1, p_2 ));

##################################################################################################################
### Generate samnple of prices at the investment horizon
N = rmvnorm(J, cbind( 0, 0 ), rbind( c(1, r), c(r, 1)));
N_1 = N[ , 1 ];
N_2 = N[ , 2 ];

U_1 = pnorm( N_1 );
U_2 = pnorm( N_2 );

aa = nu_1 / 2;
bb = 2 * s2_1;
P_1 = qgamma( U_1, aa, scale = bb);
P_2 = qlnorm( U_2, mu_2, sqrt(s2_2));

P = cbind( P_1, P_2 );

# generate sample of final wealth
W = P %*% a;

# generate sample of PnL
PnL = (P - matrix( 1, J, 1) %*% t( p )) %*% a;

# generate sample of benchmark-relative wealth
K = diag(1, 2) - p %*% t(b) / (t(b) %*% p)[1];
WRel = P %*% t(K) %*% a;

##################################################################################################################
### Plots
NumBins = round(10 * log(J));
dev.new();
plot(P_1, P_2, xlab = "P_1", ylab = "P_2" );

dev.new()
hist(W, NumBins, main = "final wealth");

dev.new();
hist(PnL, NumBins, main = "P&L");

dev.new();
hist(WRel, NumBins, main = "benchmark-relative wealth" );

