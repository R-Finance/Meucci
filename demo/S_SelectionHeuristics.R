#' Compute the r-square of selected factors, as described in A. Meucci "Risk and Asset Allocation",
#' Springer, 2005
#'
#'  @param  Who : [vector] indices for selection
#'  @param  M   : [struct] information
#'  
#'  @return g   : [scalar] r-square for the selected factors
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "SelectGoodness.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

SelectGoodness = function( Who, M )
{
	Cov_FF_k = M$Cov_FF[ Who, Who ];
	Cov_XF_k = M$Cov_XF[ , Who ];

	# abriged version of variance of error
		minCov_U = Cov_XF_k %*% (solve(Cov_FF_k) %*% matrix( Cov_XF_k ) ); 

	# abridged version of r^2
	g = sum( diag( minCov_U ) );                        

	return( g );
}

#' Naive approach for factor selection, as described in A. Meucci "Risk and Asset Allocation", Springer, 2005
#' 
#'  @param  OutOfWho : [vector] (N x 1) of selection indices
#'  @param  Metric   : [struct] metric with information on covariance
#'  
#'  @return Who      : [vector] (N x 1) indices
#'  @return Num      : [vector] (N x 1) rank of the selection 
#'  @return G        : [vector] (N x 1) r-square (cumulative)
#'
#'  @note sorted by ascending order
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "SelectNaive.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

SelectNaive = function( OutOfWho, Metric )
{
	N = length(OutOfWho);

	a = matrix( 0, 1, N );

	for( n in 1 : N )
	{
	    a[ n ] = SelectGoodness( OutOfWho[ n ], Metric );
	}

	Who = order( -a ); 

	G = matrix( NaN, N, 1);
	
	for( n in 1 : N )
	{
	    G[ n ] = SelectGoodness( Who[ 1:n ], Metric );
	}

	Num = 1 : N ;

	return( list( Who = Who, Num = Num, G = G ) )
}


#' Recursive acceptance routine for factor selection, as described in A. Meucci "Risk and Asset Allocation", Springer, 2005
#' 
#'  @param  OutOfWho : [vector] (N x 1) of selection indices
#'  @param  AcceptBy : [scalar] number of factors to accept at each iteration
#'  @param  Metric   : [struct] metric with information on covariance
#'  
#'  @return Who      : [vector] (N x 1) indices
#'  @return Num      : [vector] (N x 1) rank of the selection 
#'  @return G        : [vector] (N x 1) r-square (cumulative)
#'
#'  @note same than recursive rejection, but it starts from the empty set, instead of from the full set
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "SelectAcceptByS.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

SelectAcceptByS = function( OutOfWho, AcceptBy, Metric )
{
	N = length(OutOfWho);

	Who 	= NULL;
	Num 	= NULL;
	G   	= NULL;
	while( length(Who) < N )
	{
	    Candidates = setdiff( OutOfWho, Who );
	    
	 	if( length( Candidates )  != 1  )
	 	{
	   		Combos     = t( combn( c(Candidates), AcceptBy ) );
	    }
	    else
	    {
	    	Combos =  matrix(nchoosek(c(Candidates), AcceptBy ));
	    }

	    L = dim( Combos )[1];
	    a = matrix( 0, 1, L );
	    for( l in 1 : L )
	    {
	        a[ l ] = SelectGoodness( cbind( Who, Combos[ l, ] ), Metric );
	    }
	    g    = max( a );
	    Pick = which.max(a);
	    Who  = cbind( Who, Combos[ Pick, ] );
	    G    = cbind( G, g );
	    Num  = cbind( Num, length(Who) );
	}

	return( list( Who = Who, Num = Num, G = G ) );
}

#' Recursive rejection routine for factor selection, as described in A. Meucci "Risk and Asset Allocation", Springer, 2005
#' 
#'  @param  OutOfWho : [vector] (N x 1) of selection indices
#'  @param  RejecttBy : [scalar] number of factors to accept at each iteration
#'  @param  Metric   : [struct] metric with information on covariance
#'  
#'  @return Who      : [vector] (N x 1) indices
#'  @return Num      : [vector] (N x 1) rank of the selection 
#'  @return G        : [vector] (N x 1) r-square (cumulative)
#'
#'  @note the recursive rejection routine in Meucci (2005, section 3.4.5) to solve heuristically the above
#'     	  problem by eliminating the factors one at a time starting from the full set
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "SelectRejectByS.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

SelectRejectByS = function(OutOfWho, RejectBy, Metric)
{
    Who = OutOfWho;
    Num = length( Who );
    G   = SelectGoodness( Who, Metric );

    while( length(Who) > 1 )
    {
       
	   	Drop = t( combn( Who, RejectBy ) );
        
        L = dim( Drop )[ 1 ];
        a = matrix( 0,  1, L );
        for( l in 1 : L )
        {
            a[ l ] = SelectGoodness( setdiff( Who, Drop[ l, ] ), Metric );
        }
        g = max(a);
        Pick = which.max( a );
        Who = setdiff( Who, Drop[ Pick, ] );
        G   = cbind( G, g ); 
        Num = cbind( Num, length(Who) );
    }

    return( list( Who = Who, Num = Num, G = G ) );
}



#' Exact approach for factor selection, as described in A. Meucci "Risk and Asset Allocation", Springer, 2005
#' 
#'  @param  OutOfWho : [vector] (N x 1) of selection indices
#'  @param  Metric   : [struct] metric with information on covariance
#'  
#'  @return Who      : [vector] (N x 1) indices
#'  @return Num      : [vector] (N x 1) rank of the selection 
#'  @return G        : [vector] (N x 1) r-square (cumulative)
#'
#'  @note o iterate over the full set of factor combination
#'     	  o !!! extremely time consuming !!!
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "SelectRejectByS.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}

SelectExactNChooseK = function( OutOfWho, K, M )
{	
	Combos = t(combn( OutOfWho, K ) );
	L = dim(Combos)[1];
	a = matrix( 0, 1, L );
	
	for( l in 1 : L )
	{
	    a[ l ] = SelectGoodness( Combos[ l, ] , M );
	}

	g = max(a)
	Pick = which.max( a );
	Who  = Combos[ Pick, ];

	return( list( Who = Who, g = g ) );
}

#' This script selects the best K out of N factors in the Factors on Demand apporach to attribution
#' as described in A. Meucci, "Risk and Asset Allocation", Springer, 2005,  Chapter 3.
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170}.
#' See Meucci's script for "S_SelectionHeuristics.m"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#'  


##################################################################################################################
### Inputs
N = 50;
A = randn(N + 1, N + 1);
Sig = A %*% t(A);

Metric = list( Cov_FF = Sig[ 1:N, 1:N ], Cov_XF = matrix( Sig[ N+1, 1:N ], 1, ));
OutOfWho = 1:N;
    
##################################################################################################################
### Naive routine for factor selection
SN = SelectNaive( OutOfWho, Metric );

##################################################################################################################
### Acceptance routine for factor selection
AcceptBy = 1;
SAB = SelectAcceptByS( OutOfWho, AcceptBy, Metric );

##################################################################################################################
### Rejection routine for factor selection
RejectBy = 1;
SRB = SelectRejectByS( OutOfWho, RejectBy, Metric );

##################################################################################################################
### Plots
dev.new();
h1 = plot( SN$Num, SN$G, col = "black", type = "l", xlab = paste( "num players out of total", N ), ylab = "fit" );
h2 = lines( SAB$Num, SAB$G, col = "blue" );
h3 = lines( SRB$Num, SRB$G, col = "red" );
legend("bottomright", 1.9, c("naive", "rec. rejection", "rec. acceptance"), col = c( "black", "red", "blue" ), lty = 1, bg = "gray90" )

# exact routine
print("exact routine; be patient...");

nOutOfWho = length( OutOfWho );
GE = NULL;
NumE = NULL;

for( k in 1 : nOutOfWho )
{
    print(k);
    
    SENC = SelectExactNChooseK( OutOfWho, k, Metric );
    GE   = cbind( GE, SENC$G ); 
    NumE = cbind( NumE, k );
}

h4 = plot( NumE, GE, col = "red" );

