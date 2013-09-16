#' This script displays the sample eigenvalues dispersion phenomenon, as described in A. Meucci,
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 4.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_EigenValueDispersion.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Inputs

N = 50;
SampleLenght = seq( N , 10 * N, N)
nSim = 50;

##################################################################################################################
### Generate mesh and surface

Mu = matrix( 0, N, 1 );
Sigma= diag( 1, N );

# compute true eigenvalues
Eigen = eigen(Sigma);
Index = order( -( Eigen$values ));
EVec  = Eigen$vectors[ , Index ];
EVal  = diag( Eigen$Values[ Index, Index ]);

# compute eigenvalues of sample estimator
nSampleLenght = length( SampleLenght );
Store_EVal_Hat = matrix( NaN, nSampleLenght, N ); # preallocation for speed
for( i in 1 : nSampleLenght )
 {  
    T = SampleLenght[ i ];
    EVal_Hat = 0;
    for( n in 1 : nSim )
    {
        X = rmvnorm( T, Mu, Sigma );
        Sigma_Hat = cov( X );
        L = eigen( Sigma_Hat )$values;
        Index = order(-(L));
        L = L[ Index];

        EVal_Hat = EVal_Hat + L;
    }
    EVal_Hat = EVal_Hat / nSim;

    Store_EVal_Hat[ i, ] = t(EVal_Hat);
}

##################################################################################################################
### Display surface
dev.new();

persp( SampleLenght/N, 1 :N , Store_EVal_Hat,
    theta = 7 * 45, phi = 30, expand=0.6, col='lightblue', shade=0.75, ltheta=120, 
    ticktype='detailed', xlab = "eigenvalue #", ylab = "sample lenght/N");

