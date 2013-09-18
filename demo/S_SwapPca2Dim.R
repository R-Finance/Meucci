#' This script performs the principal component analysis of a simplified two-point swap curve.
#' it computes and plots, among others, 
#' 		1. the invariants, namely rate changes
#'
#' 		2. the location-dispersion ellipsoid of rates along with the 2-d location-dispersion ellipsoid
#'
#' 		3. the effect on the curve of the two uncorrelated principal factors 
#'
#' Described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 3. 
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 110 – Hidden factors: principal component analysis of a two-point swap curve".
#'
#' See Meucci's script for "S_SwapPca2Dim.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Load data

data("swap2y4y.mat" );

##################################################################################################################
### Current curve

Current_Curve = swap2y4y$Rates[ nrow( swap2y4y$Rates ), ];
dev.new();
plot(c( 2, 4 ), Current_Curve, type = "l", main = "Current_Curve", xlab = "time to maturity, years", ylab = "par swap rate, #" );

##################################################################################################################
### Determine weekly invariants (changes in rates)

Keep  = seq( 1, length(swap2y4y$Dates), 5 );

Rates = swap2y4y$Rates[ Keep, ];
X     = Rates[ -1, ] - Rates[ -nrow( Rates ), ];

Dates = swap2y4y$Dates[ Keep ];

PerformIidAnalysis( Dates[ -1 ], X[ , 1 ], "weekly 2yr rates" );
PerformIidAnalysis( Dates[ -1 ], X[ , 2 ], "weekly 4yr rates" );

# scatter plot of  invariants
dev.new();
plot( X[ , 1 ], X[ , 2 ], xlab = "2yr rate", ylab = "4yr rate" );
m = 0 * matrix( apply( X, 2, mean ) ); # estimator shrunk to zero
S = cov(X);
TwoDimEllipsoid(m, S, 2, TRUE, FALSE);


##################################################################################################################
### Perform PCA
E = eigen(S);
# sort eigenvalues in decreasing order
Index   = order(-E$values);
EigVals = EVals[ Index ];
EigVecs = E$vectors[ , Index ];

# plot eigenvectors
dev.new();
plot( c( 2, 4 ), EigVecs[ , 1 ], type = "l", col = "red", xlab = "Time to maturity, years", ylab = "" );
lines( c( 2, 4 ), EigVecs[ , 2 ], col = "green" );
legend("topleft", 1.9, c("1st factor loading","2nd factor loading"), col = c( "red" , "green"),
	lty = 1, bg = "gray90" );

# factors
F     = X %*% EigVecs;
F_std = apply( F, 2, sd);

dev.new(); # 1-st factor effect
plot( c( 2, 4), Current_Curve, type = "l", xlab = "Time to maturity, years", ylab = "",  ylim = c( 4.9, 5.3) );
lines( c( 2, 4), matrix( Current_Curve ) + F_std[ 1 ] * EigVecs[ , 1 ], col = "red" );
lines( c( 2, 4), matrix( Current_Curve ) - F_std[ 1 ] * EigVecs[ , 1 ], col = "green" );
legend("topleft", 1.9, c( "base", "+1 sd of 1st fact","-1 sd of 1st fact"), col = c( "black", "red" , "green"),
	lty = 1, bg = "gray90" );

dev.new(); # 2-nd factor effect
plot(  c( 2, 4), Current_Curve, type = "l", xlab = "Time to maturity, years", ylab = "",  ylim = c( 5, 5.15) );
lines( c( 2, 4), matrix( Current_Curve ) + F_std[ 2 ] * EigVecs[ , 2 ], col = "red" );
lines( c( 2, 4), matrix( Current_Curve ) - F_std[ 2 ] * EigVecs[ , 2 ], col = "green" );
legend("topleft", 1.9, c( "base", "+1 sd of 2nd fact","-1 sd of 2nd fact"), col = c( "black", "red" , "green"),
	lty = 1, bg = "gray90" );

# generalized R2
R2 = cumsum(EigVals) / sum(EigVals);  # first entry: one factor, second entry: both factors
disp(R2);

