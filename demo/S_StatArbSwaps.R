#' This script search for cointegrated stat-arb strategies among swap contracts, as described in A. Meucci, 
#' "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_StatArbSwaps.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

# TODO: Check the loadings of the principal components analysis, fix the date ticks on the plots.

##################################################################################################################
### Load data
data("swapParRates");

##################################################################################################################
### Estimate covariance and PCA decomposition
S   = cov( swapParRates$Rates );
PC  = princomp( covmat=S );
E   = PC$loadings
Lam = ( PC$sdev )^2
##################################################################################################################
### Set up dates ticks
dev.new(); 
h = plot(swapParRates$Dates, swapParRates$Dates); 
XTick = NULL;
years = as.numeric(format(swapParRates$Dates[1],"%Y")) : as.numeric(format(swapParRates$Dates[length(swapParRates$Dates)],"%Y"))


for( n in years )
{
    XTick = cbind( XTick, datenum(n,1,1) ); ##ok<AGROW>
}

a = min(swapParRates$Dates); 
b = max(swapParRates$Dates); 
X_Lim = cbind( a - 0.01 * ( b-a ),  b + 0.01 * ( b - a ) );

##################################################################################################################
### Plots
nLam = length(Lam);
Thetas = matrix( NaN, nLam, 1);
for( n in 1 : nLam )
{
    Y = swapParRates$Rates %*% E[ , n ] * 10000;
    FOU = FitOrnsteinUhlenbeck(Y, 1/252);
    Sd_Y = sqrt( FOU$Sigma / (2 * FOU$Theta));
    Thetas[n] = FOU$Theta;
    
    dev.new();
    current_line = array( Y[ length(Y) ], length(swapParRates$Dates ) );
    Mu_line = array( FOU$Mu, length(swapParRates$Dates) );
    Z_line_up = Mu_line + Sd_Y[1];
    Z_line_dn = Mu_line - Sd_Y[1];
        
    plot( swapParRates$Dates , Y, "l", xlab = "year", ylab = "basis points", main = paste( "eigendirection n. ", n, ", theta = ", FOU$Theta ) );
    lines( swapParRates$Dates, Mu_line, col = "blue" );
    lines( swapParRates$Dates, Z_line_up, col = "red" );
    lines( swapParRates$Dates, Z_line_dn, col = "red" );
    lines( swapParRates$Dates, current_line, col = "green");
    
    #set(gca(), 'xlim', X_Lim, 'XTick', XTick);
    #datetick('x','yy','keeplimits','keepticks');
    #grid off;
    #title(['eigendirection n. ' num2str(n) ',    theta = ' num2str(Theta)],'FontWeight','bold');
}

dev.new();
plot( 1:length( Lam ), Thetas, "l", xlab = " eigendirection n.", ylab = "theta" );
