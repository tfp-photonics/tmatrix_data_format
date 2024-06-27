function obj = init( obj, tau, varargin )
%  INIT - Initialize BEM solver for full Maxwell equations.
%
%  Usage for obj = galerkin.bemsolver :
%    obj = init( obj, tau, PropertyPairs )
%  Input
%    tau        :  boundary elements
%  PropertyName
%    'engine'   :  engine for potential integration

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'engine', galerkin.pot1.engine( varargin{ : } ) );
%  parse input
parse( p, varargin{ : } );

%  initialize potential integrator
obj.pot = set( p.Results.engine, tau, tau, varargin{ : } );
obj.tau = tau;
