function [ a2, b2, M2 ] = sph2pol( tab, a1, b1, t, varargin )
%  SPH2POL - Convert Mie coefficients from spherical to polar coordinates.
%
%  Usage :
%    [ a2, b2, M1 ] = sph2pol( tab, a1, b1, t, PropertyPairs )
%  Input
%    tab    :  table of spherical degrees and orders
%    a1,b1  :  Mie scattering coefficients
%    t      :  polar angle
%    ntab   :  number of tabulation points for polar angle
%  Output
%    a2,b2  :  polar scattering coefficients, farfield limit
%    M      :  angular orders

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ntab', max( 5 * max( tab.l ), 500 ) );
%  parse input
parse( p, varargin{ : } );

%  tabulation for polar angles
tt = reshape( linspace( min( t ), max( t ), p.Results.ntab ), [], 1 );
%  spherical harmonics with order M and M+1
[ L, M ] = deal( tab.l( : ), tab.m( : ) );
y2 = spharm( L, M,     tt, 0 * tt );
y1 = spharm( L, M + 1, tt, 0 * tt );

%  derivative of spherical harmonics wrt polar and azimuthal angle
y1 = M * cot( tt .' ) .* y2 + bsxfun( @times, sqrt( ( L - M ) .* ( L + M + 1 ) ), y1 );
y2 = 1i * bsxfun( @times, M, y2 );

%  prefactors
fac = ( - 1i ) .^ ( L - M + 1 ) ./ sqrt( L .* ( L + 1 ) );
%  unique angular orders and allocate output
M2 = unique( M );
[ a2, b2 ] = deal( zeros( numel( M2 ), numel( tt ), size( a1, 2 ) ) );

%  loop over unique orders
for i2 = 1 : numel( M2 )
  %  index to angular orders
  ind = M == M2( i2 ); 
  %   polar scattering coefficients, farfield limit
  a2( i2, :, : ) = bsxfun( @times, fac( ind ), y1( ind, : ) ) .' * a1( ind, : ) -  ...
               ( fac( ind ) ./ sin( tt .' ) .* y2( ind, : ) ) .' * b1( ind, : );
  b2( i2, :, : ) = bsxfun( @times, fac( ind ), y1( ind, : ) ) .' * b1( ind, : ) +  ...
               ( fac( ind ) ./ sin( tt .' ) .* y2( ind, : ) ) .' * a1( ind, : );
end
%  interpolate from table to requested polar angles
a2 = ipermute( interp1( tt, permute( a2, [ 2, 1, 3 ] ), t ), [ 2, 1, 3 ] );
b2 = ipermute( interp1( tt, permute( b2, [ 2, 1, 3 ] ), t ), [ 2, 1, 3 ] );
