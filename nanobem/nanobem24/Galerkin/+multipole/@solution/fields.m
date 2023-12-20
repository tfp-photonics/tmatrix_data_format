function [ e, h ] = fields( obj, pos, varargin )
%  FIELDS - Electromagnetic fields, Hohenester Eq. (E.4).
%
%  Usage for obj = multipole.solution :
%    [ e, h ] = fields( obj, pos )
%  Input
%    pos    :  positions where fields are computed
%  Output
%    e,h    :  electromagnetic fields

%  wavenumber and impedance
[ k, Z ] = deal( obj.mat.k( obj.k0 ), obj.mat.Z( obj.k0 ) );
%  angles of propagation directions
[ u, t, r ] = cart2sph( pos( :, 1 ), pos( :, 2 ), pos( :, 3 ) );
t = 0.5 * pi - t;
%  vector spherical harmonics
tab = obj.tab;
[ xm, xe, x0 ] = vsh( tab.l, tab.m, t, u );
%  Bessel function and derivative
x = r * k;
[ z, zp ] = riccatibessel( tab.l, x, 'h' ); 

%  dummy indices for tensor class
[ i, j, k, m ] = deal( 1, 2, 3, 4 );
%  convert positions and radial functions to tensor class
x = tensor( x, i );
z = tensor(  z,  [ m, i ] );
zp = tensor( zp, [ m, i ] );
%  prefactor
fac = tensor( sqrt( tab.l .* ( tab.l + 1 ) ), m );
xm = tensor( xm, [ m, i, k ] );
xe = tensor( xe, [ m, i, k ] );
x0 = tensor( x0, [ m, i, k ] );
%  Mie coefficients
a = tensor( obj.a, [ m, j ] );
b = tensor( obj.b, [ m, j ] );
    
%  expansion functions
M = z * xm;
N = - ( fac * z * x0 + zp * xe ) ./ x;
%  electromagnetic fields
e = ( b * M + a * N ) * Z;
h = ( a * M - b * N );
    
%  convert to numeric
e = double( sum( e, m ), [ i, k, j ] );
h = double( sum( h, m ), [ i, k, j ] );
