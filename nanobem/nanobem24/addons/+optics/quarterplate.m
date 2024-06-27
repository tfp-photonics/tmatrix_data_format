function mat = quarterplate( varargin )
%  QUARTEPLATE - Jones matrix for quarter-wavelength plate.
%
%  Usage :
%    mat = optics.quarterplate( t, PropertyPairs )
%  Input
%    t      :  angle in degress
%  PropertyName
%    dim    :  2,3 for 2d,3d matrix
%  Output
%    mat    :  Jones matrix

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 't', 0 );
addOptional( p, 'dim', 3 );
%  parse input
parse( p, varargin{ : } );

%  rotation angle
t = p.Results.t / 180 * pi;
[ sint, cost ] = deal( sin( t ), cos( t ) );
%  Jones matrix elements
fac = exp( - 0.25i * pi );
a11 = fac * ( cost ^ 2 + 1i * sint ^ 2 );
a22 = fac * ( sint ^ 2 + 1i * cost ^ 2 );
a12 = fac * ( 1 - 1i ) * sint * cost;

switch p.Results.dim
  case 2
    mat = [ a11, a12; a12, a22 ];
  case 3
    mat = [ a11, a12, 0; a12, a22, 0; 0, 0, 1 ];
end
