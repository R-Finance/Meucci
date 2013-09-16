#' This script illustrates the buy & hold dynamic strategy, as described in A. Meucci,"Risk and Asset Allocation",
#' Springer, 2005,  Chapter 6.  
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_BuyNHold.m"
#
#' @author Xavier Valls \email{flamejat@@gmail.com}

##################################################################################################################
### Input parameters
Initial_Investment = 1000;
Time_Horizon = 6 / 12; # in years
Time_Step = 1 / 252; # in years

m = 0.2;  # yearly expected return on the underlying
s = 0.40; # yearly expected percentage volatility on the stock index
r = 0.04; # risk-free (money market) interest rate

NumSimul = 30000;

##################################################################################################################
# proportion of underlying you want to hold in the beginning, e.g.: 50
Prct = 50 ;

##################################################################################################################
#### Initialize values
Underlying_Index = Initial_Investment;  # value of the underlyting at starting time, normalzed to equal investment
Start = Underlying_Index;
Elapsed_Time = 0;
Portfolio_Value = Initial_Investment;

Underlying_in_Portfolio_Percent = Prct / 100;

Underlyings_in_Portfolio = Portfolio_Value * Underlying_in_Portfolio_Percent;
Cash_in_Portfolio = Portfolio_Value - Underlyings_in_Portfolio;

##################################################################################################################
### Initialize parameters for the plot (no theory in this)

Portfolio_Series  = Portfolio_Value;
Market_Series     = Underlying_Index;
Percentage_Series = Underlying_in_Portfolio_Percent;

# asset evolution and portfolio rebalancing
while( Elapsed_Time < (Time_Horizon - 10^(-5))  ) # add this term to avoid errors
{
    # time elapses...
    Elapsed_Time = Elapsed_Time + Time_Step;
    
    # ...asset prices evolve and portfolio takes on new value...
    Multiplicator = exp( (m - s ^ 2 / 2) * Time_Step + s * sqrt( Time_Step ) * rnorm(NumSimul));
    Underlying_Index = Underlying_Index * Multiplicator;
    Underlyings_in_Portfolio = Underlyings_in_Portfolio * Multiplicator;
    Cash_in_Portfolio = Cash_in_Portfolio * exp(r * Time_Step);
    Portfolio_Value = Underlyings_in_Portfolio + Cash_in_Portfolio;
    
    # ...and we rebalance our portfolio
    Underlying_in_Portfolio_Percent = Underlyings_in_Portfolio / Portfolio_Value;    
    
    # store one path for the movie (no theory in this)
    Portfolio_Series  = cbind( Portfolio_Series, Portfolio_Value[ 1 ] ); ##ok<*AGROW>
    Market_Series     = cbind( Market_Series, Underlying_Index[ 1 ] );
    Percentage_Series = cbind( Percentage_Series, Underlying_in_Portfolio_Percent[ 1 ] );
}

##################################################################################################################
### Play the movie for one path
Time  = seq( 0, Time_Horizon, Time_Step);
y_max = max( cbind( Portfolio_Series, Market_Series) ) * 1.2;
dev.new();
par( mfrow = c(2,1))
for( i in 1 : length(Time) )
{
    plot( Time[ 1:i ], Portfolio_Series[ 1:i ], type ="l", lwd = 2.5, col = "blue", ylab = "value",
     xlim = c(0, Time_Horizon), ylim = c(0, y_max), main = "investment (blue) vs underlying (red) value");
    lines( Time[ 1:i ], Market_Series[ 1:i ], lwd = 2, col = "red" );
    #axis( 1, [0, Time_Horizon, 0, y_max]);
    
    plot(Time[ 1:i ], Percentage_Series[ 1:i ], type = "h", col = "red", xlab = "time", ylab = "#",
      xlim = c(0, Time_Horizon), ylim =c(0,1), main = "percentage of underlying in portfolio");
}

##################################################################################################################
### Plots
# plot the scatterplot
dev.new();

# marginals
NumBins = round(10 * log(NumSimul));
layout( matrix(c(1,2,2,2,1,2,2,2,1,2,2,2,0,3,3,3), 4, 4, byrow = TRUE));
barplot( table( cut( Portfolio_Value, NumBins )), horiz=TRUE, yaxt="n")

# joint scatter plot
plot(Underlying_Index, Portfolio_Value, xlab = "underlying at horizon (~ buy & hold )", ylab = "investment at horizon" );
so = sort( Underlying_Index );
lines( so, so, col = "red" );

barplot( table( cut( Underlying_Index, NumBins )), yaxt="n")
