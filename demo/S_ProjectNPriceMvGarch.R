#'This script fits a multivariate GARCH model and projects the distribution of the compounded returns 
#' from the estimation interval to the investment horizon. 
#'Then it computes the distribution of prices at the investment horizon , as described in A. Meucci,
#'"Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_ProjectNPriceMvGarch.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

if ( !require( "signal" ) ) stop("signal package installation required for this script")

##################################################################################################################
### Load data
data("equities" );

##################################################################################################################
### Inputs

Prices = Equities$Prices[ , c(4, 5)];
J = 10000; # numbers of MC scenarios
N = ncol(Prices); # numbers of securities 
T = 22; # projection horizon 

##################################################################################################################
### Estimation of daily compounded returns distribution
CR = diff(log(Prices)); # extract risk drivers (compounded returns)
demean = 1;
eps    = .01;
df     = 500;
#[m, A, B, C, H] 
MVGarch = FitMultivariateGarch(CR, demean, eps, df);

##################################################################################################################
### Projection to monthly compouded returns distribution
H_ = array( 0, dim = c( J, N, N ))
for( j in 1 : J )
{
    H_[ j, ,  ] = MVGarch$H;
}

X_T = matrix( 0, J, N );
for( t in 1 : T )
{
    for( j in 1 : J )# WARNING: this loop is for didactical purposes only. In real applications avoid looping
    {
        # compute new return
        e = matrix(rnorm(N));
        H = drop( H_[ j, , ] );
        X = MVGarch$mu + chol(H) %*% e;
        X_T[ j, ] = X_T[ j, ] + t(X);
    
        # update for next cycle
        S = X %*% t(X);
        H = MVGarch$CTMF + MVGarch$ATMF * S + MVGarch$BTMF * H;
        H_[ j, , ] = H;
    }
}

##################################################################################################################
### Pricing into linear returns distribution
R = exp( X_T ) - 1;

##################################################################################################################
### Display results
dev.new()
# marginals
NumBins = round( 10 * log( J ) );

layout( matrix(c(1,2,2,1,2,2,0,3,3), 3, 3, byrow = TRUE), heights=c(1,2,1));
# marginal
barplot( table( cut( X_T[ , 2 ], NumBins )), horiz=TRUE, yaxt="n")

# scatter plot
plot( X_T[ , 1 ], X_T[ , 2 ], xlab = "", ylab = "" );

# marginal
barplot( table( cut( X_T[ , 1 ], NumBins )), yaxt="n")
