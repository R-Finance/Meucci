#' Interpolate and extrapolate using n-linear interpolation (tensor product linear).
#'
#'  @param   V        : [array] p-dimensional array to be interpolated/extrapolated at the list of points in the array Xi.
#                       interpne will work in any number of dimensions >= 1
#'  @param   Xi       : [array] (n x p) array of n points to interpolate/extrapolate. Each point is one row of the array Xi.
#'  @param   nodelist : [cell array] (optional) cell array of nodes in each dimension.
#                       If nodelist is not provided, then by default I will assume nodelist[[i]] = 1:size(V,i). The nodes in
#                       nodelist need not be uniformly spaced.
#'  @param   method   : [string] (optional) chacter string, denotes the interpolation method used. default method = 'linear'
#                       'linear'  --> n-d linear tensor product interpolation/extrapolation
#                       'nearest' --> n-d nearest neighbor interpolation/extrapolation
#                       in 2-d, 'linear' is equivalent to a bilinear interpolant
#                       in 3-d, it is commonly known as trilinear interpolation.
#'  
#'  @return Vpred     : [array] (n x 1) array of interpolated/extrapolated values
#'  
#'  @note   
#'  Initially written by John D'Errico
#'  Vpred = interpne(V,Xi)
#'  Vpred = interpne(V,Xi,nodelist)
#'  Vpred = interpne(V,Xi,nodelist,method)
#'  Extrapolating long distances outside the support of V is rarely advisable.
#'
#' @references
#' \url{http://symmys.com/node/170}
#' See Meucci's script for "InterExtrapolate.R"
#'
#' @author Xavier Valls \email{flamejat@@gmail.com}
#' @export

#  examples
#
#  [x1,x2] = meshgrid(0:.2:1);
#  z = exp(x1+x2);
#  Xi = rand(100,2)*2-.5;
#  Zi = interpne(z,Xi,{0:.2:1, 0:.2:1},'linear');
#  surf(0:.2:1,0:.2:1,z)
#  hold on
#  plot3(Xi(:,1),Xi(:,2),Zi,'ro')
#

InterExtrapolate = function( V, Xi, nodelist, method )
{
    # get some sizes

    vsize = dim( V );
    ndims = length( vsize );
    nargin = length( match.call() ) -1;
    if( ndims != ncol( Xi ) ) stop("Xi is not compatible in size with the array V for interpolation.")

    # default for nodelist
    if ( nargin < 3 || length(nodelist) == 0 )
    {
        nodelist= vector( "list", ndims );

        for( i in 1 : ndims )
        {
            nodelist[[ i ]] = rbind( 1 : vsize[ i ] );
        }
    }

    if ( length( nodelist ) != ndims ) 
        stop( "nodelist is incompatible with the size of V.")

    nll = lapply( nodelist, length );

    if( any( nll!= vsize) ) 
        stop( "nodelist is incompatible with the size of V." )

    # get deltax for the node spacing
    dx = nodelist;

    for( i in 1 : ndims )
    {
        nodelist[[i]] = nodelist[[i]][]; # not sure about this doing anything
        dx[[i]] = diff(nodelist[[i]][,]);
        if ( any( dx[[i]] <= 0) )
            stop( "The nodes in nodelist must be monotone increasing." );
    }

    # check for method
    if ( nargin < 4 ) method = 'linear';

    

    if( ! is.character(method) ) stop("method must be a character string if supplied.");
   
    
    validmethod = c( "linear", "nearest");
    
    if(!any(validmethod == method ) ) 
        stop(paste(" No match found for method =  ", method))

    # Which cell of the array does each point lie in?
    # This includes extrapolated points, which are also taken
    # to fall in a cell. histc will do all the real work.

    ind = matrix( 0, nrow(Xi), ndims);
    
    for( i in 1 : ndims)
    {
        hc = histc(Xi[ , i ], nodelist[[i]]); ##ok<ASGLU>
         
        # catch any point along the very top edge.
        hc$bin[ hc$bin == vsize[ i ] ] = vsize[ i ] - 1;
        
        ind[ , i ] = hc$bin;
        
        k = which( hc$bin == 0);
           
        # look for any points external to the nodes
        if( !(length(k)==0) )
        {
            # bottom end
            ind[ k[ Xi[ k, i] < nodelist[[i]][ 1 ]], i ] = 1;
            
            # top end
            ind[ k[ Xi[ k, i] > nodelist[[i]][ length(nodelist[[1]][1]) ]], i ] = vsize[ i ] - 1;
        }
    }

    # where in each cell does each point fall?
    t = matrix( 0, nrow(Xi), ndims);

    for( i in 1 : ndims)
    {
        t[ , i ] = (Xi[ , i ] - nodelist[[i]][ ind[ , i ] ] )/ dx[[i]][ind[ , i ]];
    }

    sub = cumprod( c( 1 ,vsize[ 1 : ( length( vsize ) - 1 ) ] ) );
    base = 1 + ( ind-1 ) %*% sub;

    # which interpolation method do we use?

    switch( method,
        nearest = 
        {
            # nearest neighbor is really simple to do.
            t = round(t);
            t[ t > 1 ] = 1;
            t[ t < 0 ] = 0;
            
            Vpred = V[ base + t %*% sub ];
        },

        linear = 
        {
            # tensor product linear is not too nasty.
            Vpred = matrix( 0, nrow(Xi), 1);
            # define the 2^ndims corners of a hypercube (MATLAB's corners = (dec2bin(0:(2^ndims-1))== '1');)
            corners = lapply( strsplit( intToBin ( 0 : ( 2^ndims - 1 ) ), split=""), as.integer ); 
            
            nc = length( corners );
            
            for( i in 1 : nc )
            {
                #accessing
                s = V[ base + (corners[[i]] %*% sub)[1]]; 
                for( j in 1 : ndims )
                {
                    # this will work for extrapolation too
                    if( corners[[i]][ j ] == 0 ){
                        s = s * ( 1 - t[ , j ] );
                    }else
                    {
                        s = s * t[ , j ];
                    }
                }

                Vpred = Vpred + s;
            }
        }   ) # end switch method
    
    return( Vpred );
}