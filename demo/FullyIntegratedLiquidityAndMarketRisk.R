#' This script computes the liquidity-risk and funding-risk adjusted P&L distribution, as described in 
#' A. Meucci, "A Fully Integrated Liquidity and Market Risk Model", Financial Analyst Journal, 68, 6, 35-47 (2012)
#'
#' @references 
#' \url{http://www.symmys.com/node/350}
#' See Meucci script "S_Main.m"
#' 
#' @author Xavier Valls \email{flamejat@@gmail.com}

# INPUTS
#####################################################################*
# liquidation policy at horizon as fraction of investment
Policy=-1;

# collinearity of liquidity perturbations
CollinLiq=1;

# select only some stock in portfolio and equally allocate capital as fraction of daily dollar volume
Selectstock = 1:10 ;
Capital_perDailyVolume = 0.2;


# PREPARE DATA
#####################################################################*
# load  fILMR$Daily_Prices: closing prices 
#       fILMR$Daily_Volumes_Shares: daily volumes   
#       fILMR$Daily_Liq: Morgan Stanley liquidity index 
data("fILMR")

# Prices and returns
#Daily_Prices = Daily_Prices(:,Selectstock);
Prices_0 = matrix(  fILMR$Daily_Prices[ nrow(fILMR$Daily_Prices), ] );
Daily_LogRets = log( fILMR$Daily_Prices[ -nrow(fILMR$Daily_Prices), ] / fILMR$Daily_Prices[ -1, ] );
J = dim( Daily_LogRets )[1]
N = dim( Daily_LogRets )[2]

# volumes in shares
#Daily_Volumes = Daily_Volumes_Shares[ , Selectstock ];
Volumes_0 =  matrix( fILMR$Daily_Volumes_Shares[ nrow(fILMR$Daily_Volumes_Shares), ]);

Volumes_t = matrix( apply( fILMR$Daily_Volumes_Shares[ -(1:(nrow(fILMR$Daily_Volumes_Shares)-250)), ], 2, mean ) );

# liquidity index
Daily_LiqChanges = diff(fILMR$Daily_Liq);
Liq_0 = matrix( fILMR$Daily_Liq[ length(fILMR$Daily_Liq), ] );

# normal simulations
    X    = cbind( Daily_LogRets, Daily_LiqChanges );
    m_X  = apply( X, 2, mean );
    s2_X = cov( X ); #covariance
    J    = 100000;
    #RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', 11)); 
    X = rmvnorm( J, m_X, s2_X );  
    Daily_LogRets    = X[ ,1:N ];
    Daily_LiqChanges = X[ , dim(X)[2] ];

# Fully Flexible Probabilties associated with each scenario
Probs = matrix( 1, J, 1 ) / J;

# stock prices at horizon 
Prices_t = repmat( t( Prices_0) , J, 1 ) * exp(Daily_LogRets);

# liquidity index at horizon
Liq_t = Liq_0 * exp(Daily_LiqChanges);

# pure market risk: p&L due to market risk
PnL_mkt = Prices_t - repmat( t(Prices_0), J, 1 );


# PORTFOLIO COMPUTATIONS
######################################################################
# portfolio and liquidation policy
Weights = matrix( 0, N, 1);
Weights[ Selectstock ] = 1 / length(Selectstock);
DollarVolume_0 = t(Volumes_0) %*% Prices_0;
Capital = (Capital_perDailyVolume %*% DollarVolume_0)[1];

h = Capital * Weights / Prices_0;

PnL_mkt_h = PnL_mkt %*% h;                           

# LIQUIDITY ADJUSTMENT 
######################################################################
# liquidation policy
Dh = Policy * h;

# market impact
b_a       = 0.01 * matrix( 1, N, 1 );                          
Linear    =-b_a * Prices_0 * abs(Dh);
NonLinear = -(10^5) * Prices_0 * matrix( apply( Daily_LogRets, 2, sd ) ) * ( ( abs(Dh)  / Volumes_t ) ^ 1.5);         
m_Dh      = Linear + NonLinear;                                      
    
# state-dependent expected liquidity impact on all stocks
s_g1   = 0.5 * sd( PnL_mkt_h );
g1     = -pmin( PnL_mkt_h, -s_g1 ) / s_g1;
m_Dh_x = repmat( g1, 1, N ) * repmat( t(m_Dh), J, 1 );    # (14)

# state-dependent expected liquidity impact on portfolio
m_Dh_h = m_Dh_x %*% matrix( 1, N, 1 );                    # (23)

# state-independent uncertainty on liquidity impact on portfolio
s_Dh    = 1.5 * m_Dh;                                                                       # 
r2_Dh   = ( 1 - CollinLiq ) * cor( Daily_LogRets ) + CollinLiq * matrix( 1, N, N );         # 
s2_Dh   = diag( s_Dh[,] , length(s_Dh) ) %*% r2_Dh %*% diag( s_Dh[,], length( s_Dh ) );     # 
s2_Dh_h = t( matrix( 1, N, 1 ) ) %*% s2_Dh %*% matrix( 1, N, 1 );                           # 
s_Dh_h  = sqrt(s2_Dh_h);    
s_Dh_h  = pmax( s_Dh_h, 0.01 * std(PnL_mkt_h) );                                            # regularization

# TOTAL P&L
######################################################################
# conditional center and scatter
m_j = PnL_mkt_h + m_Dh_h;
s_j = s_Dh_h[1] * matrix( 1, J, 1 );

# pdf and cdf: taking and not taking into account funding cost
nu   = 100;
f_Pi = function(x){ t( Probs / s_j) %*% dt( (x-m_j) / s_j, nu ) };
F_Pi = function(x){ t( Probs ) %*% pt( (x-m_j) / s_j, nu ) };

NGrid = 200;
x_= seq( min(PnL_mkt_h) - s_Dh_h, max(PnL_mkt_h) + s_Dh_h, length = NGrid );
p_= NULL;
f_Pi_plot = NULL;
f_Pi_funding_plot = NULL;
for( k in 1:NGrid )
{
    p_= rbind( p_, F_Pi( x_[ k ] ) );
    f_Pi_plot = rbind( f_Pi_plot, f_Pi( x_[ k ] ) );
}

###########################################################################
# plots
dev.new()
NumBins = round( 10 * log(J) );
hist = hist(PnL_mkt_h,NumBins, plot = FALSE ); # compute bin width
D = hist$mids[ 2 ]- hist$mids[ 1 ];
hh = plot( hist$mids, hist$counts / ( J * D ), type = "h", xlab ="", ylab = "" ); # plot histogram
hh1 = lines(x_,f_Pi_plot, col="red");

legend( "topright", 1.9, c("pure market P&L","market + liquidity P&L"), col = c( "black", "red"), lty=1, bg = "gray90" );