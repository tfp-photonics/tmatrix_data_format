function [ p, poly ] = tripolygon( poly, edge, varargin )
%  TRIPOLYGON - Three-dimensional particle from polygon.
%
%  Usage :
%    [ p, poly ] = tripolygon( poly, edge, PropertyPairs )
%  Input
%    poly        :  2d polygon
%    edge        :  edge profile
%  PropertyName
%    'hdata'     :  pass hdata to mesh2d
%    'opt'       :  pass options structure to mesh2d
%    'refine'    :  refine function for edges
%  Output
%    p3         :  particle object for extruded particle
%    poly       :  edge polygon

%  rounded or sharp upper and lower edges
if all( ~isnan( edge.pos( :, 1 ) ) ) || ~( nnz( edge.pos( :, 1 ) == 0 ) == 1 )
  %  polygon3 object
  poly1 = arrayfun( @( poly ) polygon3( poly, edge.zmin ), poly, 'UniformOutput', false );
  poly2 = arrayfun( @( poly ) polygon3( poly, edge.zmax ), poly, 'UniformOutput', false );
  %  plates
  [ plate1, ~    ] = plate( horzcat( poly1{ : } ), 'edge', edge, 'dir', -1, varargin{ : } );
  [ plate2, poly ] = plate( horzcat( poly2{ : } ), 'edge', edge, 'dir',  1, varargin{ : } );
  %  ribbon
  ribbon = vribbon( poly, varargin{ : } );
  
%  sharp lower edge  
elseif isnan( edge.pos( 1, 1 ) )
  %  polygon3 object
  poly = arrayfun( @( poly ) polygon3( poly, edge.zmax ), poly, 'UniformOutput', false );  
  %  upper plate
  [ plate1, poly1 ] = plate( horzcat( poly{ : } ), 'edge', edge, 'dir',  1, varargin{ : } );
  %  ribbon
  [ ribbon, ~, poly ] = vribbon( poly1, varargin{ : } );  
  %  lower plate
  plate2 = fun( plate1, poly1, set( poly, 'z', edge.zmin ) );
  
%  sharp upper edge  
else
  %  polygon3 object
  poly = arrayfun( @( poly ) polygon3( poly, edge.zmin ), poly, 'UniformOutput', false );  
  %  lower plate
  [ plate1, poly1 ] = plate( horzcat( poly{ : } ), 'edge', edge, 'dir', -1, varargin{ : } );
  %  ribbon
  [ ribbon, poly, ~ ] = vribbon( poly1, varargin{ : } );  
  %  upper plate
  plate2 = fun( plate1, poly1, set( poly, 'z', edge.zmax ) );
  plate2.faces = fliplr( plate2.faces );
  
end

%  put together particle
ribbon = num2cell( ribbon );
p = union( plate1, plate2, ribbon{ : } ); 



function plate2 = fun( plate1, poly1, poly2 )
%  FUN - Shift plate circumference from POLY1 to POLY2.

%  polygon positions
pos1 = vertcat( poly1.poly.pos );
pos2 = vertcat( poly2.poly.pos );
%  find vertices at circumference
[ verts, faces ] = deal( plate1.verts( :, 1 : 2 ), plate1.faces );
[ i1, i2 ] = ismember( verts, pos1, 'rows' );
%  shift plate circumference and smooth mesh
verts( i1, : ) = pos2( i2( i1 ), : );
[ verts, faces ] = smoothmesh( verts, faces, 1000 );
%  set particle
verts( :, 3 ) = poly2.z;
plate2 = particle( verts, faces );
