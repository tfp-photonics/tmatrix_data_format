function y = cumeval( obj, r, z, fun, k0, varargin )
%  CUMEVAL - Cumulative Sommerfeld integral (for testing only).
%
%  Usage for obj = stratified.isommerfeld :
%    y = cumeval( obj, r, z, fun, k0, PropertyPairs )
%  Input
%    r,z        :  radii and z-values for evaluation 
%    fun        :  Sommerfeld integrand
%    k0         :  wavenumber of light in vacuum
%  PropertyName
%    singular   :  singular value subtraction of Chew
%    nquad      :  number of integration points
%  Output
%    y          :  cumulative Sommerfeld integral

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'singular', 0 );
addParameter( p, 'nquad', 100 );
%  parse input
parse( p, varargin{ : } );

switch p.Results.singular
  case 0
    y = cumeval1( obj, r, z, fun, k0, p.Results.nquad );
  case 1
    y = cumeval2( obj, r, z, fun, k0, p.Results.nquad );
end
