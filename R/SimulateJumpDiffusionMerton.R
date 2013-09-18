#' @title Simulates a Merton jump-diffusion process.
#'
#' @description This function simulates a jump diffusion process, as described in A. Meucci "Risk and Asset Allocation",
#' Springer, 2005.
#'
#'  @param  m   [scalar] deterministic drift of diffusion
#'  @param  s   [scalar] standard deviation of diffusion
#'  @param  l   [scalar] Poisson process arrival rate
#'  @param  a   [scalar] drift of log-jump
#'  @param  D   [scalar] st.dev of log-jump
#'  @param  ts  [vector] time steps
#'  @param  J   [scalar] number of simulations
#'
#'  @return  X  [matrix] (J x length(ts)) of simulations
#'
#' @references
#' A. Meucci - "Exercises in Advanced Risk and Portfolio Management" \url{http://symmys.com/node/170},
#' "E 132 - Simulation of a jump-diffusion process".
#'
#' See Meucci's script for "SimulateJumpDiffusionMerton.m"
#'
#' Merton, R. C., 1976. "Option pricing when underlying stocks are discontinuous". Journal of Financial
#' Economics 3, 125-144.
#' 
#'@author Xavier Valls \email{flamejat@@gmail.com}
#' @export

SimulateJumpDiffusionMerton = function( m, s, l, a, D, ts, J )
{
    L = length(ts);
    T = ts[ L ];

    # simulate number of jumps; 
    N = rpois( J, l * T );

    Jumps = matrix( 0, J, L );
    for( j in 1 : J )
    {
        # simulate jump arrival time
        t = T * rnorm(N[ j ]);
        t = sort(t);

        # simulate jump size
        S = a + D * rnorm(N[ j ]);
        
        # put things together
        CumS = cumsum(S);
        Jumps_ts = matrix( 0, 1, L);
        for( n in 1 : L )
        {
            Events = sum( t <= ts[ n ]);
            if( Events )
            {
                Jumps_ts[ n ] = CumS[ Events ];
            }
        }

        Jumps[ j, ] = Jumps_ts;
    }

    D_Diff = matrix( NaN, J, L );
    for( l in 1 : L )
    {
        Dt = ts[ l ];
        if( l > 1 )
        {
            Dt = ts[ l ] - ts[ l - 1 ];
        }

        D_Diff[ , l ] = m * Dt + s * sqrt(Dt) * rnorm(J); 
    }

    X = cbind( matrix(0, J, 1), apply(D_Diff, 2, cumsum) + Jumps );

    return( X );
}