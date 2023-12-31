function obj = init( obj, varargin )
%  INIT - Initialize default evaluation of quasistatic BEM potentials.
%
%  Usage for obj = galerkinstat.pot1.std :
%    obj = init( obj, PropertyPairs )
%  PropertyName
%    'rules'      :  quadrature rules

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules', quadboundary.rules );
%  parse input
parse( p, varargin{ : } );

obj.parent = true;
obj.rules = p.Results.rules;
