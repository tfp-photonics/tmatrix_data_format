function y = eval( obj, r, z, fun, k0, varargin )
%  EVAL - Evaluate Sommerfeld integral.
%
%  Usage for obj = stratified.isommerfeld :
%    y = eval( obj, r, z, fun, k0, PropertyPairs )
%  Input
%    r,z        :  radii and z-values for evaluation 
%    fun        :  Sommerfeld integrand
%    k0         :  wavenumber of light in vacuum
%  PropertyName
%    singular   :  singular value subtraction of Chew
%  Output
%    y          :  Sommerfeld integral

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'singular', 0 );
%  parse input
parse( p, varargin{ : } );

switch p.Results.singular
  case 0
    y = eval1( obj, r, z, fun, k0 );
  case 1
    y = eval2( obj, r, z, fun, k0 );
end
