function obj = init( obj, layer, pol, dir, varargin )
%  INIT - Initialize plane wave object for layer structure.
%
%  Usage for obj = galerkin.planewave :
%    obj = init( obj, pol, varargin )
%  Input
%    pol    :  polarizations of plane waves

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules', quadboundary.rules );
%  parse input
parse( p, varargin{ : } );

%  layer structure, plane wave polarization and light propagation direction
obj.layer = layer;
obj.pol = pol;
obj.dir = dir;
%  integration rules 
obj.rules = p.Results.rules;
