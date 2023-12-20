function obj = init1( obj, layer, i1, varargin )
%  INIT1 - Initialize base class for tabulated Green functions.
%
%  Usage for obj = stratified.tab.base :
%    obj = init1( obj, layer, i1, PropertyPairs )
%  Input
%    layer    :  layer structure
%    i1       :  medium index
%  PropertyPairs 
%    method   :  method for interpolation
      
%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'method', 'cubic' );
%  parse input
parse( p, varargin{ : } );

%  layer structure and interpolation method
obj.layer = layer;
obj.method = p.Results.method;
%  initialize Sommerfeld integrator
obj.somint = stratified.isommerfeld( layer, i1, varargin{ : } );
