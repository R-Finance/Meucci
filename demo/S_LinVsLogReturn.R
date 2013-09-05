#' This script project a distribution in the future according to the i.i.d.-implied square-root rule, as described
#'  in A. Meucci "Risk and Asset Allocation", Springer, 2005, chapter 3.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "S_LinVsLogReturn.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}


##################################################################################################################
### Inputs
# in general R=exp(C)-1. Furtheremore, here we assume C~N(m*t,s^2*t) 
m = 0.05;
s = 0.25;

ts = seq( 0.1, 3, 0.3 );
ps = matrix( c( 0.01, 0.99 ) );

D = 0.7 * min( abs( diff( ts ) ) );
C = list( q = NULL, x = NULL, pdf = NULL);
R = list( q = NULL, x = NULL, pdf = NULL);

Steps = 100;

for( i in 1 : length(ts) )
{
    t = ts[ i ];
    
    q = qnorm( ps, m * t ,s * sqrt( t ) );
    
    x = seq( min(q), max(q) , (max(q)-min(q))/Steps );

    pdf = dnorm( x, m*t, s*sqrt(t));
    pdf = pdf / max( pdf ) * D;
    
    C$q   = cbind( C$q, q );
    C$pdf = cbind( C$pdf, pdf );
    C$x   = cbind( C$x, x );

    q = exp( q )-1;
    
    x = seq( min(q), max(q), (max(q)-min(q))/Steps );

    pdf = dlnorm( x + 1, m * t, s * sqrt( t ) );
    pdf = pdf / max( pdf ) * D;
    
    R$pdf = cbind( R$pdf, pdf);
    R$x   = cbind( R$x, x );
}
 

R$q = exp( C$q ) - 1;

Col = rgb( 0.8, 0.8, 0.8 );


par(mfrow=c(2,1));


matplot(c( 0, ts ), t(cbind( 0*ps, C$q )), type="l", lty=1, col = "red", 
    xlab ="", ylab ="", main = "compounded returns" );


for( i in 1 : length(ts) )
{
    xx = rbind( ts[i] , ts[i] + C$pdf[ ,i ] , ts[i]);
    yy = rbind( min(C$x[,i]) , C$x[ ,i ], max(C$x[,i]));
    polygon(xx, yy, col= Col);
}


matplot(c( 0, ts ), t(cbind( 0*ps, R$q )), type="l", lty=1, col = "red", 
    xlab ="", ylab ="", main = "linear returns" );

for( i in 1 : length(ts) )
{
    xx = rbind( ts[i] , ts[i] + R$pdf[ ,i ] , ts[i]);
    yy = rbind( min(R$x[,i]) , R$x[ ,i ], max(R$x[,i]));
    polygon(xx, yy, col= Col);
}

# xlim and ylim  should be ylim = min(yy)*1.5 max(yy)*1.5,  xlim=c(0, max(xx)*1.01)) respectively in each one of the plots