function [ e, h ] = farfields( obj, dir, varargin )
%  FARFIELDS - Scattered farfields, Hohenester Eq. (E.27).
%
%  Usage for obj = multipole.solution :
%    [ e, h ] = farfields( obj, dir, PropertyPairs )
%  Input
%    dir    :  propagation direction
%  PropertyName
%    shift  :  additional shift for far-fields
%  Output
%    e,h    :  electromagnetic farfields

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'shift', [] );
%  parse input
parse( p, varargin{ : } );

%  wavenumber and impedance
[ k, Z ] = deal( obj.mat.k( obj.k0 ), obj.mat.Z( obj.k0 ) );
%  angles of propagation directions
[ u, t ] = cart2sph( dir( :, 1 ), dir( :, 2 ), dir( :, 3 ) );
t = 0.5 * pi - t;
%  vector spherical harmonics
tab = obj.tab;
[ xm, xe ] = vsh( tab.l, tab.m, t, u );

%  multiply scattering coefficients with prefactor
a = bsxfun( @times, obj.a, ( - 1i ) .^ ( tab.l( : ) + 1 ) );
b = bsxfun( @times, obj.b, ( - 1i ) .^ ( tab.l( : ) + 1 ) );
%  multiplication function
n1 = size( xm, 1 );
n2 = size( xm, 2 );
fun = @( a, x ) permute(  ...
  reshape( a .' * reshape( x, n1, [] ), [], n2, 3 ), [ 2, 3, 1 ] );
%  electromagnetic farfields, Hohenester Eq. (E.27)
e = squeeze( fun( b, xm ) - 1i * fun( a, xe ) ) / k * Z;
h = squeeze( fun( a, xm ) + 1i * fun( b, xe ) ) / k;

%  deal with shift argument
if ~isempty( p.Results.shift )
  %  phase factor
  fac = exp( - 1i * k * dir * p.Results.shift .' );
  %  multiply electromagnetic fields with phase factor
  e = double( tensor( e, [ 1, 2, 3 ] ) * tensor( fac, [ 1, 3 ] ), [ 1, 2, 3 ] );
  h = double( tensor( h, [ 1, 2, 3 ] ) * tensor( fac, [ 1, 3 ] ), [ 1, 2, 3 ] );
end
