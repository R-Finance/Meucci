#' This script shows that the eigenvectors of a Toeplitz matrix have a Fourier basis structure under t-distribution
#' assumptions, as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 130 – Eigenvectors for Toeplitz structure".
#'
#' See Meucci's script for "S_Toeplitz.R"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}


###############################################################################################################
### Inputs

N = 200; # dimension of the matrix
Decay = 0.9; # decay factor

###############################################################################################################
T = diag( 1, N);
for( n in 1 : (N - 1) )
{

    T = T + Decay^n * ( cbind(  matrix( 0, N, N -( N-n ) ), diag( 1, N , N-n) ) +
      cbind( rbind( matrix(0, N-(N-n), N-n ), diag( 1, N-n)), matrix(0, N, N-(N-n) ) )) ;

}
eig = eigen( T );

###############################################################################################################

#R sorts the eigen vectors, so the results aren't going to be exactly the same as in MATLAB

dev.new();
plot( eig$vectors[ , n ], type = "l", col = runif(1)*100 );
lines( eig$vectors[ , n-1 ], type = "l", col = runif(1)*100 );
