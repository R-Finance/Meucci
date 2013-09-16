
#'This script studies horizon effect on explicit factors / implicit loadings linear model, as described in 
#'A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'Compounded returns follow the linear model X = tau*muX + D*F + epsilon, where 
#' tau: investment horizon (in weeks)
#' muX: expected weekly compounded returns
#' F: factor compounded returns, with zero expectation and tau-proportional covariance
#' D: matrix of factor loadings
#' epsilon: uncorrelated (idiosyncratic) shocks.
#' R = exp(X)-1 and Z = exp(F)-1 are the linear returns
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_HorizonEffect.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
# Load parameters of the model: D, muX, sigmaF, sigmaEps
data("linearModel" );

# Specify range of investment horizon, weeks
tauRangeWeeks = 1:52;

# constants

N = nrow(linearModel$D);
K = ncol(linearModel$D);

ntauRangeWeeks = length(tauRangeWeeks);
aMinusTauMuX = matrix(0, ntauRangeWeeks);
normDminusB  = matrix(0, ntauRangeWeeks);
minCorrU     = matrix(0, ntauRangeWeeks);
meanCorrU    = matrix(0, ntauRangeWeeks);
maxCorrU     = matrix(0, ntauRangeWeeks);

for( i in 1 : ntauRangeWeeks )
{
    tau = tauRangeWeeks[ i ];
    
    # statitics of linear returns (analytically)
    # let Y = [1+R; 1+Z] ~ LogN(mu*tau, covJointXF*tau)
    mu = rbind(linearModel$muX, matrix(0, K));
    # covariance of [X; F] for tau=1:
    covJointXF = rbind(
        cbind( linearModel$D %*% linearModel$sigmaF %*% t(linearModel$D) + linearModel$sigmaEps,
         linearModel$D %*% linearModel$sigmaF ),
        cbind( linearModel$sigmaF %*% t( linearModel$D ), linearModel$sigmaF)
        );
    E_Y  = exp( mu * tau + diag( covJointXF * tau ) / 2 );
    E_YY = (E_Y %*% t( E_Y )) * exp( covJointXF * tau );
    E_R  = E_Y[ 1 : N ] - matrix( 1, N );
    E_Z  = E_Y[ -(1:N)] - matrix( 1, K);
    E_RR = E_YY[ 1:N , 1:N ] - matrix( 1, N ) %*% t(E_R) - E_R %*% matrix( 1, 1, N ) - matrix( 1, N, N );
    E_ZZ = E_YY[-(1:N), -(1:N) ] - matrix( 1, K ) %*% t(E_Z) - E_Z %*% matrix( 1, 1, K ) - matrix( 1, K, K );
    E_RZ = E_YY[ 1:N, -(1:N) ] - matrix( 1, N ) %*% t(E_Z) - E_R %*% matrix( 1, 1, K) - matrix( 1, N, K );
    SigmaZ  = E_ZZ - E_Z %*% t(E_Z);
    SigmaR  = E_RR - E_R %*% t(E_R);
    SigmaRZ = E_RZ - E_R %*% t(E_Z);

    # compute OLS loadings for the linear return model
    B = SigmaRZ %*% solve( SigmaZ ); # right division of SigmaRZ by SigmaZ
    a = E_R - B %*% E_Z;
    aMinusTauMuX[ i] = norm( a - tau*linearModel$muX, type="2" );
    normDminusB[ i ]  = norm( linearModel$D - B, type = "F" );
    
    # pairwise correlations of U
    SigmaU = SigmaR - B %*% t( SigmaRZ );
    corrU = cov2cor(SigmaU);
    stackedCorrU = as.array(corrU);
    minCorrU[ i ]  = min(min(abs(corrU)));
    meanCorrU[ i ] = ( N * mean( mean( abs( corrU) ) ) - 1 ) / ( N - 1 );
    corrU[ seq(1, N * N, (N+1) ) ] = 0;
    maxCorrU[ i ] = max( max( abs( corrU) ) );
}

##################################################################################################################
### Plots
# relationship between the constant nd the intercept
dev.new();
plot(tauRangeWeeks, aMinusTauMuX, type= "l", xlab = expression(paste("investment horizon, ", tau,", weeks")),
 main = expression( paste( "norm of ( a - ", tau,mu[X],")"^t )));

# relationship between the loadings D in and the loadings B
dev.new();
plot(tauRangeWeeks, normDminusB, type = "l", xlab = expression(paste("investment horizon, ", tau,", weeks")), main = expression("norm of (D-B)"^t));



# determine if U idiosyncratic
dev.new();
plot(tauRangeWeeks, maxCorrU, col = "red", type = "l", xlab = expression(paste("investment horizon, ", tau,", weeks")),
    ylab = "", main = "pairwise correlations of elements of U" );
lines(tauRangeWeeks, meanCorrU, col = "blue");
lines(tauRangeWeeks, minCorrU, col = "green");
legend( "topleft", 1.9, c( "max absolute corr", "mean absolute corr", "min absolute corr" ), col = c( "red","blue", "green" ),
     lty=1, bg = "gray90" );
