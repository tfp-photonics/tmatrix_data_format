function [ e, h, data ] = fields( obj, pt1, varargin )
%  FIELDS - Electromagnetic fields at point positions.
%
%  Usage for obj = stratified.solution :
%    [ e, h, data ] = fields( obj, pt1, PropertyPairs )
%  Input
%    pt1    :  point positions
%  PropertyName
%    pot1   :  direct potential integrator
%    pot2   :  reflected potential integrator
%    green  :  reflected Green function object
%    refl   :  compute reflected fields only
%  Output
%    e      :  electric field at requested points
%    h      :  magnetic field at requested points
%    data   :  potential integrators and Green function for reuse

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'pot1',  [] );
addParameter( p, 'pot2',  [] );
addParameter( p, 'green', [] );
addParameter( p, 'refl',   0 );
%  parse input
parse( p, varargin{ : } );

%  extract input
[ pot1, pot2, green ] =  ...
  deal( p.Results.pot1, p.Results.pot2, p.Results.green );

%  observation points and boundary elements directly connected ?
ind = bsxfun( @eq,  ...
  reshape( vertcat( obj.tau.inout ), [], 1 ), horzcat( pt1.imat ) );
%  direct potential integrator
if isempty( pot1 ) && any( ind( : ) )
  pot1 = galerkin.pot2.engine( varargin{ : } );
  pot1 = set( pot1, pt1, obj.tau, varargin{ : } );
end
%  reflected potential integrator
if isempty( pot2 )
  pot2 = stratified.pot2.reflected( obj.layer, pt1, obj.tau, varargin{ : } );
end
%  reflected Green function
if isempty( green )
  pos2 = arrayfun( @( q ) eval( q ), quadboundary( obj.tau ), 'uniform', 0 );
  pos2 = reshape( cat( 1, pos2{ :  } ), [], 3 );
  r = slice( obj.layer, vertcat( pt1.pos ), pos2 );
  r = grid( obj.layer, r, varargin{ : } );
  green = stratified.tab.green( r, varargin{ : } );
end

%  fill Green function
green = fill( green, obj.k0 );
%  direct fields, consider case that boundary elements and observation
%  points are not directly connected
n = numel( pt1 );
if isempty( pot1 ) || p.Results.refl
  [ e1, h1 ] = deal( 0 );
else
  [ e1, h1 ] = fields( pot1, obj, varargin{ : }, 'n', n );
end
%  reflected fields
[ e2, h2 ] = fields( pot2, green, obj, 'n', n );
%  set output
e = e1 + e2;
h = h1 + h2;
data = struct( 'pot1', pot1, 'pot2', pot2, 'green', green );
