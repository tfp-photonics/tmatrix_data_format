function obj = init( obj, tau, layer, varargin )
%  INIT - Initialize BEM solver for layer structures.
%
%  Usage for stratified.bemsolver :
%    obj = init( obj, tau, layer, PropertyPairs )
%  Input
%    tau        :  boundary elements
%    layer      :  layer structure
%  PropertyName
%    'engine'   :  engine for direct potential integrator

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'engine', galerkin.pot1.engine( varargin{ : } ) );
%  parse input
parse( p, varargin{ : } );

%  save boundary elements and layer structure
[ obj.tau, obj.layer ] = deal( tau, layer );
%  reflected Green function table
r = slice( layer, vertcat( tau.pos ), vertcat( tau.pos ) );
r = grid( layer, r, varargin{ : } );
obj.green = stratified.tab.green( r, varargin{ : } );
%  initialize potential integrators
obj.pot1 = set( p.Results.engine, tau, tau, varargin{ : } );
obj.pot2 = stratified.pot1.reflected( layer, tau, tau, varargin{ : } );
