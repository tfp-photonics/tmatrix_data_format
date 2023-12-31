function [ dmin, imin, pmin ] = dist( obj, pos )
%  DIST - Find distance of point positions from polygon.
%
%  Usage for obj = polygon :
%    [ dmin, imin, pmin ] = dist( obj, pos )
%  Input
%    obj    :  single object or array
%    pos    :  positions(s)
%  Output
%    dmin   :  minimum distance of each position to polygons
%    imin   :  nearest polygon vertex for each position
%    pmin   :  pointer to nearest polygon object

dmin =  ones( size( pos, 1 ), 1 ) * 1e10;
imin = zeros( size( pos, 1 ), 1 );
pmin = zeros( size( pos, 1 ), 1 );

%  loop over polygons
for i = 1 : numel( obj )

  %  positions of polygon
  xa = obj( i ).pos( :, 1 );  xb = obj( i ).pos( [ 2 : end, 1 ], 1 );
  ya = obj( i ).pos( :, 2 );  yb = obj( i ).pos( [ 2 : end, 1 ], 2 );

  %  loop over polygon positions
  for it = 1 : size( pos, 1 )
    
    x = pos( it, 1 );
    y = pos( it, 2 );
  
    lambda = ( ( xb - xa ) .* ( x  - xa ) + ( yb - ya ) .* ( y  - ya ) ) ./  ...
             ( ( xb - xa ) .* ( xb - xa ) + ( yb - ya ) .* ( yb - ya ) );
    lambda = max( min( lambda, 1 ), 0 );
  
    d = sqrt( ( xa + lambda .* ( xb - xa ) - x ) .^ 2 +  ...
              ( ya + lambda .* ( yb - ya ) - y ) .^ 2 );
    [ d, ind ] = min( d );
    
    if d < dmin( it )
      dmin( it ) = d;        %  minimum distance
      imin( it ) = ind;      %  nearest position      
      pmin( it ) = i;        %  nearest polygon
    end
  
  end
end
