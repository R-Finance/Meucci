#' This script dispays location-dispersion ellipsoid and statistic, as described 
#' in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_MaxMinVariance.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

if ( !require( "mvtnorm" ) ) stop("mvtnorm package installation required for this script")

##################################################################################################################
### Input parameters
Mu   = rbind( 0.5, 0.5 );
s    = rbind( 0.1, 0.1 );
nu   = 40;
r    = -0.9;
nSim = 10000;

##################################################################################################################
### Generate sample
C    = rbind( c( 1, r ), c( r, 1));
Ones = matrix( 1, nSim, 1);
Y    = Ones %*% t(Mu) + ( Ones %*% t(s) ) * rmvt( nSim, C, nu );
X    = exp( Y );
m    = matrix( apply( X, 2, mean ));
S    = cov( X );

##################################################################################################################
### Evaluate standard deviation on a one-dim projection (versor)
Theta  = seq( 0, 2 * pi, pi/100 );
nTheta = length( Theta );
invS = solve( S );

SDev   = matrix( NaN, nTheta, 1 );
Radius = matrix( NaN, nTheta, 1 );

for( n in 1 : nTheta )
{
    th = Theta[ n ];
    # versor
    e = rbind( cos(th), sin(th) );
    Z = X %*% e;  # projection
    SDev[ n ]   = sd( Z );  # standard deviation    
    Radius[ n ] = ( t(e) %*% invS %*% e )^( -1/2 );   # radius of ellipsoid
}

##################################################################################################################
### Compute min and max standard deviation and respective versor
min_n = which.min( SDev );
e_min = rbind( cos( Theta[ min_n ] ), sin( Theta[ min_n ] ) );
s_min = SDev[ min_n ];

max_n = which.max( SDev );
e_max = rbind( cos( Theta[ max_n ] ), sin( Theta[ max_n ] ) );
s_max = SDev[ max_n ];

##################################################################################################################
### Plots
dev.new();

# scatter plot simulations
plot( X[ , 1 ], X[ , 2 ] );

# plot ellipsoid
Scale = 2;
PlotEigVectors = 1;
PlotSquare = 1;
TwoDimEllipsoid( m, S, Scale, PlotEigVectors, PlotSquare);

# plot special directions defined by the max-min versors
Center = matrix( apply( X, 2, mean) ) * 0.7; # de-center plot of special directions for better display
v = Scale * matrix( seq( -1, 1, 0.1 ) );
Ones = 1 + 0 * v;

v_min = Ones %*% t( Center ) + v %*% s_min %*% t( e_min );
lines( v_min[ , 1 ], v_min[ , 2 ], col = "red" );
v_max = Ones %*% t( Center ) + v %*% s_max %*% t( e_max );
lines( v_max[ , 1 ], v_max[ , 2 ], col = "red" );

# plot statistics versus geometry
dev.new();
Scaled_Theta = Theta / (pi / 2);
 # plot standard deviation as function of direction
plot( Scaled_Theta, SDev, type = "l", xlab = "theta/(pi/2)", xlim = c( Scaled_Theta[ 1 ], Scaled_Theta[length(Scaled_Theta)] ) );
# plot radius of ellipsoid as function of direction
lines( Scaled_Theta, Radius, col="red" ); 
legend( "topleft", 1.9, c( "st.dev on projection", "radius of ellipsoid" ), col = c( "black", "red" ), lty = 1, bg = "gray90" );