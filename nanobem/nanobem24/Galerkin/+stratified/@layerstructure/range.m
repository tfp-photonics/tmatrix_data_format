function r = range( obj, pts, varargin )
%  RANGE - Tabulation ranges for position pairs.
%
%  Usage for obj = stratified.layerstructure :
%    r = range( obj, pts, PropertyPairs )
%  Input
%    pts      :  position pairs, see SLICE
%  PropertyName
%    rmin     :  minimum radial distance
%    zmin     :  minimum distance to interface
%    margin   :  enlarging factor for bounding box
%  Output
%    r        :  range structure

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rmin', 1e-2 );
addParameter( p, 'margin', 0.05 );
%  parse input
parse( p, varargin{ : } );

%  deal with input vector
if ~isscalar( pts )
  r = arrayfun( @( x ) range( obj, x, varargin{ : } ), pts, 'uniform', 1 );
else
  %  layer
  r.layer = obj;
  %  layer indices
  r.i1 = pts.i1;
  r.i2 = pts.i2;
  
  %  maximum and minimum distance between positions
  dmin = pdist2( pts.pos1( :, 1 : 2 ),  ...
                 pts.pos2( :, 1 : 2 ), 'euclidean', 'smallest', 1 );
  dmax = pdist2( pts.pos1( :, 1 : 2 ),  ...
                 pts.pos2( :, 1 : 2 ), 'euclidean', 'largest',  1 );
  %  minimum and maximum of radial distance
  if isscalar( dmin )
    r.r = bsxfun( @plus, p.Results.rmin * [ -1, 1 ], dmin );
  else
    r.r = [ max( min( dmin ), p.Results.rmin ), ( 1 + p.Results.margin ) * max( dmax ) ];  
  end
   
  %  z-range
  z1 = [ min( pts.pos1( :, 3 ) ), max( pts.pos1( :, 3 ) ) ];
  z2 = [ min( pts.pos2( :, 3 ) ), max( pts.pos2( :, 3 ) ) ];
  if z1( 1 ) == z1( 2 ),  z1( 2 ) = z1( 2 ) + 1e-2;  end
  if z2( 1 ) == z2( 2 ),  z2( 2 ) = z2( 2 ) + 1e-2;  end
  %  observation and source points in same layer
  if r.i1 == r.i2
    switch r.i1
      case 1
        %  lowest layer
        zlim = bsxfun( @minus, 2 * obj.z( 1 ), [ z1( 2 ) + z2( 2 ),  ...
                                                 z1( 1 ) + z2( 1 ) ] );
        r.z1 = fun( zlim, [ 0, 1e30 ], varargin{ : } );
        r.z2 = [];            
      case obj.n + 1
        %  highest layer
        zlim = bsxfun( @minus, [ z1( 1 ) + z2( 1 ),  ...
                                 z1( 2 ) + z2( 2 ) ], 2 * obj.z( end ) );
        r.z1 = fun( zlim, [ 0, 1e30 ], varargin{ : } );
        r.z2 = [];   
      otherwise
        %  inside layer structure
        z1 = fun( z1, obj.z( r.i1 + [ -1, 0 ] ), varargin{ : } );
        z2 = fun( z2, obj.z( r.i1 + [ -1, 0 ] ), varargin{ : } );
        %  sum and difference of z-values
        [ Z1, Z2 ] = deal( obj.z( r.i1 - 1 ), obj.z( r.i1 ) );
        r.z1 = [ z1( 1 ) - z2( 2 ), z1( 2 ) - z2( 1 ) ] + Z2 - Z1;
        r.z2 = [ z1( 1 ) + z2( 1 ), z1( 2 ) + z2( 2 ) ] - 2 * Z1;
    end
  else
    z = [ - 1e30, reshape( obj.z, 1, [] ), 1e30 ];
    r.z1 = fun( z1, z( r.i1 + [ 0, 1 ] ), varargin{ : } );
    r.z2 = fun( z2, z( r.i2 + [ 0, 1 ] ), varargin{ : } );    
  end
end


function z = fun( z, zlayer, varargin )
%  FUN - Set z-range.

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'zmin', 1e-2 );
addParameter( p, 'margin', 0.1 );
%  parse input
parse( p, varargin{ : } );

dz = p.Results.margin * diff( z );
z = [ max( z( 1 ) - dz, zlayer( 1 ) + p.Results.zmin ),   ...
      min( z( 2 ) + dz, zlayer( 2 ) - p.Results.zmin ) ];  
