function obj = init( obj, tau, param, varargin )
%  INIT - Initialize BEM with Feibelman parameters.
%
%  Usage :
%    obj = feibelman.bemsolver( tau, param, PropertyPairs )
%  Input
%    tau    :  boundary elements
%    param  :  Feibelman parameters
%  PropertyName
%    'engine'   :  potential integrator

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'engine', galerkin.pot1.engine( varargin{ : } ) );
%  parse input
parse( p, varargin{ : } );

%  initialize potential integrator
obj.pot = set( p.Results.engine, tau, tau, varargin{ : } );
obj.tau = tau;
%  indices for connected boundary elements, Feibelman parameters
obj.ind = connectivity( tau );
obj.param = param;
