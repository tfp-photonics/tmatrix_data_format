function [ d, iz ] = bdist2( obj, varargin )
%  BDIST2 - Pairwise distance to layer scaled by bounding box radius.
%
%  Usage for obj = stratified.layerstructure :
%    [ d, iz ] = bdist2( obj, tau1, tau2 )
%    [ d, iz ] = bdist2( obj, tau1, pos2 )
%  Input
%    tau    :  boundary elements
%    pos    :  positions
%  Output
%    d      :  distance scaled by bounding box radius
%    iz     :  layer index for closest interface

switch class( varargin{ 2 } )
  case 'BoundaryEdge'
    %  extract input
    [ tau1, tau2 ] = deal( varargin{ : } );
    %  centroids of boundary elements and bounding boxes
    pos1 = vertcat( tau1.pos );  rad1 = boxradius( tau1 ); 
    pos2 = vertcat( tau2.pos );  rad2 = boxradius( tau2 );  
    %  minimal distance to interface
    [ d, iz ] = fun( obj, pos1, pos2 );
    d = d ./ bsxfun( @plus, rad1, rad2 .' );
  otherwise
     %  extract input
     [ tau1, pos2 ] = deal( varargin{ : } );    
     %  centroids of boundary elements and bounding boxes
     pos1 = vertcat( tau1.pos );  rad1 = boxradius( tau1 ); 
     %  minimal distance to interface
     [ d, iz ] = fun( obj, pos1, pos2 );
     d = bsxfun( @rdivide, d, rad1( : ) );
end



function [ d, iz ] = fun( obj, pos1, pos2 )
%  FUN - Minimal distance between POS1 and POS2.

%  layer indices
ind1 = indlayer( obj, pos1 );  
ind2 = indlayer( obj, pos2 );  
%  distance to given interface of layer
dist = @( z, ind ) abs( z - obj.z( ind ) );
%  number of interfaces and allocate output
n = obj.n;
[ d, iz ] = deal( nan( numel( ind1 ), numel( ind2 ) ) );

%  loop over unique layer indices
for i1 = unique( reshape( ind1, 1, [] ) )
for i2 = unique( reshape( ind2, 1, [] ) )
  %  radial distance
  k1 = ind1 == i1;
  k2 = ind2 == i2;
  r = pdist2( pos1( k1, 1 : 2), pos2( k2, 1 : 2 ) );
  %  z-values
  z1 = pos1( k1, 3 );
  z2 = pos2( k2, 3 ) .';
  %  find minimal distance to interfaces
  if i1 == 1 && i2 == 1
    d( k1, k2 ) = r + bsxfun( @plus, dist( z1, 1 ), dist( z2, 1 ) );
    iz( k1, k2 ) = 1;
  elseif i1 == n + 1 && i2 == n + 1
    d( k1, k2 ) = r + bsxfun( @plus, dist( z1, n ), dist( z2, n ) );
    iz( k1, k2 ) = n;
  elseif i1 == i2
    d1 = bsxfun( @plus, dist( z1, i1 - 1 ), dist( z2, i1 - 1 ) );
    d2 = bsxfun( @plus, dist( z1, i1     ), dist( z2, i1     ) ); 
    d( k1, k2 ) = r + bsxfun( @min, d1, d2 );
    iz( k1, k2 ) =  ...
      ( i1 - 1 ) * bsxfun( @le, d1, d2 ) + i1 * bsxfun( @gt, d1, d2 );
  elseif i1 == i2 + 1
    d( k1, k2 ) = r + bsxfun( @plus, dist( z1, i2 ), dist( z2, i2 ) );
    iz( k1, k2 ) = i1;
  elseif i2 == i1 + 1
    d( k1, k2 ) = r + bsxfun( @plus, dist( z1, i1 ), dist( z2, i1 ) );
    iz( k1, k2 ) = i2;
  end
end
end
