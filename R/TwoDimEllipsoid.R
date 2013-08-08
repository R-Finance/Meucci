#' This script computes the location-dispersion ellipsoid of the normalized (unit variance, zero expectation)
#' first diagonal and off-diagonal elements of a 2x2 Wishart distribution as a function of the inputs,
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 2.
#'
#'  @param	Location 	      : [vector] (2 x 1) location vector (typically the expected value
#'	@param	Square_Dispersion : [matrix] (2 x 2) scatter matrix Square_Dispersion (typically the covariance matrix)
#'  @param	Scale             : [scalar] a scalar Scale, that specifies the scale (radius) of the ellipsoid
#'  @param	PlotEigVectors    : [boolean] true then the eigenvectors (=principal axes) are plotted
#'  @param	PlotSquare        : [boolean] true then the enshrouding box is plotted. If Square_Dispersion is the covariance
#'
#'	@return	E                 : [figure handle]
#'
#' @references
#' \url{http://}
#' See Meucci's script for "TwoDimEllipsoid.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export


TwoDimEllipsoid = function( Location, Square_Dispersion, Scale = 1, PlotEigVectors = FALSE, PlotSquare = FALSE )
{

	##########################################################################################################
	### compute the ellipsoid in the r plane, solution to  ((R-Location)' * Dispersion^-1 * (R-Location) ) = Scale^2                                   
	
	Eigen = eigen(Square_Dispersion);
	Centered_Ellipse = c(); 
	Angle = seq( 0, 2*pi, pi/500 );
	NumSteps = length(Angle);
	
	for( i in 1 : NumSteps )
	{
	    # normalized variables (parametric representation of the ellipsoid)
	    y = rbind( cos( Angle[ i ] ), sin( Angle[ i ] ) );
	    Centered_Ellipse = c( Centered_Ellipse, Eigen$vectors %*% diag(sqrt(Eigen$values), length(Eigen$values)) %*% y );   ##ok<AGROW>
	}

	R = Location %*% array( 1, NumSteps ) + Scale * Centered_Ellipse;

	##########################################################################################################
	### Plot the ellipsoid
	
	E = lines( R[1, ], R[2, ], col = "red", lwd = 2 );

	##########################################################################################################
	### Plot a rectangle centered in Location with semisides of lengths Dispersion[ 1]  and Dispersion[ 2 ], respectively
	
	if( PlotSquare )
	{
	    Dispersion = sqrt( diag( Square_Dispersion ) );
	    Vertex_LowRight_A = Location[ 1 ] + Scale * Dispersion[ 1 ]; 
	    Vertex_LowRight_B = Location[ 2 ] - Scale * Dispersion[ 2 ];
	    Vertex_LowLeft_A  = Location[ 1 ] - Scale * Dispersion[ 1 ]; 
	    Vertex_LowLeft_B  = Location[ 2 ] - Scale * Dispersion[ 2 ];
	    Vertex_UpRight_A  = Location[ 1 ] + Scale * Dispersion[ 1 ]; 
	    Vertex_UpRight_B  = Location[ 2 ] + Scale * Dispersion[ 2 ];
	    Vertex_UpLeft_A   = Location[ 1 ] - Scale * Dispersion[ 1 ]; 
	    Vertex_UpLeft_B   = Location[ 2 ] + Scale * Dispersion[ 2 ];
	    
	    Square = rbind( c( Vertex_LowRight_A, Vertex_LowRight_B ), 
	              c( Vertex_LowLeft_A,  Vertex_LowLeft_B ),
	              c( Vertex_UpLeft_A,   Vertex_UpLeft_B ),
	              c( Vertex_UpRight_A,  Vertex_UpRight_B ),
	              c( Vertex_LowRight_A, Vertex_LowRight_B ) );

	        h = lines(Square[ , 1 ], Square[ , 2 ], col = "red", lwd = 2 );
	        
	}

	##########################################################################################################
	### Plot eigenvectors in the r plane (centered in Location) of length the square root of the eigenvalues (rescaled)
	if( PlotEigVectors )
	{
	    L_1 = Scale * sqrt( Eigen$values[ 1 ] );
	    L_2 = Scale * sqrt( Eigen$values[ 2 ] );
	    
	    # deal with reflection: matlab chooses the wrong one
	    Sign = sign( Eigen$vectors[ 1, 1 ] );

	    # eigenvector 1
	    Start_A = Location[ 1 ];                               
	    End_A   = Location[ 1 ] + Sign * (Eigen$vectors[ 1, 1 ]) * L_1;
	    Start_B = Location[ 2 ];
	    End_B   = Location[ 2 ] + Sign * (Eigen$vectors[ 1, 2 ]) * L_1;
	    
	    h = lines( c( Start_A, End_A ), c( Start_B, End_B ), col = "red", lwd = 2 );
	    
	    # eigenvector 2
	    Start_A = Location[ 1 ];                               
	    End_A   = Location[ 1 ] + ( Eigen$vectors[ 2, 1 ] * L_2);
	    Start_B = Location[ 2 ];
	    End_B   = Location[ 2 ] + ( Eigen$vectors[ 2, 2 ] * L_2);
	    
	    h = lines( c( Start_A, End_A ), c( Start_B, End_B ), col = "red", lwd = 2 );
	    
	}
}