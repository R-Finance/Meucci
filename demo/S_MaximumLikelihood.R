#' This script performs ML under a non-standard parametric set of distributions, as described in A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_MaximumLikelihood.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

## Parametric pdf used in the ML estimation
fparam = function( x,  th )
{
    

    m = th;
        
    if( th <= 0 )
    {
        s = sqrt(th^2);
        nu = 1;
        f = 1 / s * dt((x-m) / s, nu);
    } else
    {
        s = sqrt((th - 0.01) ^ 2);
        f = dlnorm(x, m, s);
    }

    return( f );
}


qparam = function( p, th )
{
    ## Parametric inverse cdf used in the ML estimation
    m = th;
        
    if( th <= 0 )
    {
        s = sqrt(th^2);
        nu = 1;
        q = m + s * qt(p, nu);
    } else
    {
        s = sqrt((th - 0.01) ^ 2);
        q = qlnorm(p, m, s);
    }

    return( q );
}


##########################################################################################################
### Load data
data("timeSeries");

##########################################################################################################
### inputs
p = 0.01;
Theta = cbind( seq( -0.04, -0.01, 0.001 ), 0.02, 0.03 );

##########################################################################################################
### check invariance
T = length(TimeSeries$i_T);
PerformIidAnalysis( 1:T, TimeSeries$i_T, "Invariance Analysis" );

##########################################################################################################
### ML-estimate parameters
nTheta = length(Theta);
Store_LL = matrix( NaN, nTheta, 1); # preallocation for speed
for( s in 1 : nTheta )
{
    print(s);
    theta = Theta[ s ];
    # compute log-likelihood
    LL = 0;
    for( t in 1 : T )
    {
        x_t = TimeSeries$i_T[ t ];
        LL = LL + log( fparam( x_t, theta ) );
    }
    Store_LL[ s ] = LL;
}

# compute log-likelihood faster (if likelihood function is vectorized)
Store_LL_fast = matrix( NaN, nTheta, 1); # preallocation for speed
for( s in 1 : nTheta )
{
    print(s);  
    Store_LL_fast[ s ] = sum( log ( fparam( TimeSeries$i_T, Theta[ s ] ) ) );
}

# comparison
print( cbind( Store_LL, Store_LL_fast ) );

# determine the maximum likelihood
Max_Index = which.max( Store_LL );
# and the corresponding estimator
theta_ML = Theta[ Max_Index ];
print(theta_ML);

# display the LL value for range of parameters
dev.new();
plot( Theta, Store_LL, type = "o");

# compute MLE-implied quantile
Q_ML = qparam( p, theta_ML );

# compute sample quantile
Q_NP = quantile( TimeSeries$i_T, p );

# comparison of quantiles
print(cbind( Q_ML, Q_NP ));

