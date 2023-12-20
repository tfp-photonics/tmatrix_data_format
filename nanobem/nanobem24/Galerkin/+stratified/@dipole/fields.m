function [ e, h ] = fields( obj, pt, k0, varargin )
%  FIELDS - Electromagnetic fields for dipole.
%
%  Usage for obj = stratified.dipole :
%    [ e, h ] = fields( obj, pt, k0, PropertyPairs )
%  Input
%    pt     :  positions where field is evaluated
%    k0     :  wavelength of light in medium
%  PropertyName
%    green  :  tabulated Green function
%    refl   :  reflected fields only
%  Output
%    e      :  electric field
%    h      :  magnetic field

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'green', [] );
addParameter( p, 'refl', 0 );
%  parse input
parse( p, varargin{ : } );

%  field and dipole positions and material indices
pos1 = vertcat(     pt.pos );  ind1 = vertcat(     pt.imat );
pos2 = vertcat( obj.pt.pos );  ind2 = vertcat( obj.pt.imat );
%  reflected Green function
green = p.Results.green;
if isempty( green )
  r = slice( obj.layer, pos1, pos2 );
  r = grid( obj.layer, r, obj.opt{ : } );
  green = stratified.tab.green( r, obj.opt{ : } );
end
%  dyadic Green function object for reflected Green function
green = fill( green, k0 );
green = stratified.dyadic( green );

%  direct contribution
switch p.Results.refl
  case 0
    [ e, h ] = fields( obj.dip, pt, k0 );
  case 1
    [ e, h ] = deal( zeros( numel( pt ), 3, numel( obj.pt ), 3 ) );
end
%  reflected contribution
[ i1, i2 ] = deal( ind1 <= obj.layer.n + 1, ind2 <= obj.layer.n + 1 );
data = eval( green, pos1( i1, : ), pos2( i2, : ), k0, varargin{ : } );
%  assign output
e( i1, :, i2, : ) = e( i1, :, i2, : ) +  k0 ^ 2 * data.SL1;
h( i1, :, i2, : ) = h( i1, :, i2, : ) + 1i * k0 * data.DL2;
