function [ M, N ] = transvec( obj, pos, k0, varargin )
%  TRANSVEC - Transverse vector functions, Hohenester Eq. (D.3).
%
%  Usage for obj = multipole.base :
%    [ M, N ] = transvec( obj, pos, k0, PropertyPairs )
%  Input
%    pos    :  positions
%    k0     :  wavenumber
%  PropertyName
%    conj   :  complex conjugate of VSHs
%  Output
%    M,N    :  transverse vector function

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'conj', 0 );
%  parse input
parse( p, varargin{ : } );

%  convert positions to spherical coordinates
[ u, t, r ] = cart2sph( pos( :, 1 ), pos( :, 2 ), pos( :, 3 ) );
t = 0.5 * pi - t;
%  spherical degrees and orders
[ ltab, mtab ] = deal( obj.tab.l, obj.tab.m );
%  Bessel function and vector spherical harmonics
x = r * k0;
[ z, zp ] = riccatibessel( ltab, x, 'j' ); 
[ xm, xe, x0 ] = vsh( ltab, mtab, t, u );

%  dummy indices for tensor class
[ i, k, m ] = deal( 1, 2, 3 );
%  convert arrays to tensor class
[ z, zp ] = deal( tensor( z, [ m, i ] ), tensor( zp, [ m, i ] ) );
[ xm, xe, x0 ] = deal( tensor( xm, [ m, i, k ] ),  ...
                       tensor( xe, [ m, i, k ] ), tensor( x0, [ m, i, k ] ) );
%  prefactor
fac = tensor( sqrt( ltab .* ( ltab + 1 ) ), m );
x = tensor( x, i );
%  expansion functions
M = z * xm;
N = - ( fac * z * x0 + zp * xe ) ./ x;

%  set output
M = double( M, [ m, i, k ] );
N = double( N, [ m, i, k ] );
%  complex conjugate
if p.Results.conj
  [ M, N ] = deal( conj( M ), conj( N ) );
end
