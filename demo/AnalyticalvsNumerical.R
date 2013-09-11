#' This example script compares the numerical and the analytical solution of entropy-pooling, as described
#' in A. Meucci, "Fully Flexible Views: Theory and Practice", The Risk Magazine, October 2008, p 100-106.
#'
#' Most recent version of article and MATLAB code available at
#' http://www.symmys.com/node/158
#'
#' @references 
#' A. Meucci, "Fully Flexible Views: Theory and Practice" \url{http://www.symmys.com/node/158}
#' See Meucci script for "AnalyticalvsNumerical/S_MAIN.m"
#' 
#' @author Ram Ahluwalia \email{ram@@wingedfootcapital.com} and Xavier Valls \email{flamejat@@gmail.com}

###############################################################
# prior
###############################################################

# analytical representation

N  = 2 # market dimension (2 assets)
Mu = matrix( 0, N , 1 )
r  = 0.6
Sigma = ( 1 - r ) * diag( 1, N ) + r * matrix( 1, N , N ) # nxn correlation matrix with correlation 'r' in off-diagonals

# numerical representation
J  = 100000 # number of scenarios
p  = matrix( 1, J , 1 ) / J
dd = rmvnorm( J / 2 , matrix( 1, N , 1 ) , Sigma )      # distribution centered on (0,0) with variance Sigma
X  = matrix( 1, J , 1 ) %*% t(Mu) + rbind( dd , -dd )   # JxN matrix of scenarios

###############################################################
# views
###############################################################

# location
Q = cbind(  1 , -1  ) # long the first and asset and short the second asset produces an expectation (of Mu_Q calculated below)
Mu_Q = 0.5

# scatter
G = cbind(  -1 , 1  )
Sigma_G = 0.5^2

###############################################################
#  posterior 
###############################################################

# analytical posterior
RevisedMuSigma = Prior2Posterior( Mu , Q , Mu_Q , Sigma , G , Sigma_G )
Mu_ = RevisedMuSigma$M_
Sigma_ = RevisedMuSigma$S_
# numerical posterior
Aeq = matrix( 1, 1 , J )  # constrain probabilities to sum to one...
beq = 1

# create views

QX = X %*% t(Q) # a Jx1 matrix

Aeq = rbind( Aeq , t(QX) )    # ...constrain the first moments... 
    # QX is a linear combination of vector Q and the scenarios X

beq = rbind( beq , Mu_Q )

SecMom = G %*% Mu_ %*% t(Mu_) %*% t(G) + Sigma_G  # ...constrain the second moments... 
    # We use Mu_ from analytical result. We do not use Revised Sigma because we are testing whether
    # the numerical approach for handling expectations of covariance matches the analytical approach
    # TODO: Can we perform this procedure without relying on Mu_ from the analytical result?

GX = X %*% t(G)

for ( k in 1:nrow( G ) )
{
    for ( l in k:nrow( G ) )
    {
        Aeq = rbind( Aeq , t(GX[ , k ] * GX[ , l ] ) )
        beq = rbind( beq , SecMom[ k , l ] )
    }
}

emptyMatrix = matrix( , nrow = 0 , ncol = 0 )
p_ = EntropyProg( p , emptyMatrix , emptyMatrix , Aeq , beq ) # ...compute posterior probabilities

###############################################################
# plots
###############################################################
PlotDistributions( X , p , Mu , Sigma , p_ , Mu_ , Sigma_ )