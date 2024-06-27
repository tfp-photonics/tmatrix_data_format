function [ e, h ] = fields2( obj, pts, layer, varargin )
%  FIELDS2 - Fields for planewave decomposition and layer structure.
%
%  Usage for obj = optics.decompose :
%    [ e, h ] = fields2( obj, pts, layer, PropertyPairs )
%  Input
%    pts      :  points where electromagnetic fields are requested
%    layer    :  layer structure
%  PropertyName
%    primary  :  reflected fields only
%  Output
%    e,h      :  electromagnetic fields

%  index to points connected to layer structure
ind = vertcat( pts.imat ) <= layer.n + 1;
%  allocate output
[ e, h ] = deal( zeros( numel( pts ), 3 ) );

if nnz( ind )
  %  electric and magnetic fields
  pos = vertcat( pts( ind ).pos );
  [ e1, h1 ] = planewave( layer, obj.efield, obj.dir, obj.k0, pos, varargin{ : } );
  %  update fields
  e( ind, : ) = sum( e1, 3 );
  h( ind, : ) = sum( h1, 3 );
end
