function [ e, h ] = farfields( obj, pt, k0 )
%  FARFIELDS - Electromagnetic far-fields for dipole.
%
%  Usage for obj = galerkin.dipole :
%    [ e, h ] = farfields( obj, pt, k0 )
%  Input
%    pt     :  positions where field is evaluated
%    k0     :  wavelength of light in medium
%  Output
%    e      :  electric far-field
%    h      :  magnetic far-field

%  material properties and indices of points
mat = pt( 1 ).mat;
imat = horzcat( pt.imat );
%  wavenumber, refractive index, and impedance in media
k1 = arrayfun( @( x ) x.k( k0 ), mat, 'uniform', 1 );  k1 = k1( imat );
n1 = arrayfun( @( x ) x.n( k0 ), mat, 'uniform', 1 );  n1 = n1( imat );
Z1 = arrayfun( @( x ) x.Z( k0 ), mat, 'uniform', 1 );  Z1 = Z1( imat );

%  dummy indices for internal tensor class
[ i1, i2, j, k ] = deal( 1, 2, 3, 4 );
%  dipole positions, orientations and propagation directions
pos = tensor( vertcat( obj.pt.pos ), [ i2, k ] );
dip = tensor( eye( 3 ), [ j, k ] );
dir = tensor( vertcat( pt.pos ), [ i1, k ] );
%  prefactor
[ k1, n1, Z1 ] = deal( tensor( k1, i1 ), tensor( n1, i1 ), tensor( Z1, i1 ) );
fac = exp( - 1i * k1 * dot( pos, dir, k ) );

%  electromagnetic far-fields, Jackson Eq. (9.20)
h = fac * k1 .^ 2 ./ ( 4 * pi * n1 ) * cross( dir, dip, k );
e = cross( h, dir, k ) * Z1;
%  points and dipoles connected ?
is = tensor( connected( pt, obj.pt ), [ i1, i2 ] );
e = e * is;
h = h * is;
%  convert to numeric
e = double( e( i1, k, i2, j ) );
h = double( h( i1, k, i2, j ) );
