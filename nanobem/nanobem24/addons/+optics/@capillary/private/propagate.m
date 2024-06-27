function [ a, b ] = propagate( obj, a, b, k0, kz )
%  PROPAGATE - Propagate polar coefficients through capillary.
%
%  Usage for obj = optics.capillary :
%    [ a, b ] = propagate( obj, a, b, kz )
%  Input
%    a,b    :  polar coefficients
%    k0     :  wavenumber of light in vacuum
%    kz     :  z-components of wavevectors
%  Output
%    a,b    :  polar coefficients in outside medium
%    M      :  angular orders

%  reshape coefficients
siz = size( a );
a = reshape( permute( a, [ 2, 1, 3 ] ), numel( kz ), [] );
b = reshape( permute( b, [ 2, 1, 3 ] ), numel( kz ), [] );
%  material properties and cylinder diameters
[ mat, diameter ] = deal( obj.mat, obj.diameter );
%  Fresnel coefficients
[   ~, t12 ] = fresnel( mat( 1 ), mat( 2 ), k0, kz );
[ r21, ~   ] = fresnel( mat( 2 ), mat( 1 ), k0, kz );
[ r23, t23 ] = fresnel( mat( 2 ), mat( 3 ), k0, kz );
%  modulus and radial component of wavevectors
k1 = mat( 1 ).k( k0 );  kr1 = sqrt( k1 ^ 2 - kz .^ 2 );
k2 = mat( 2 ).k( k0 );  kr2 = sqrt( k2 ^ 2 - kz .^ 2 );
k3 = mat( 3 ).k( k0 );  kr3 = sqrt( k3 ^ 2 - kz .^ 2 );
%  propagation constants
psi2 = 0.5 * kr2 * ( diameter( 2 ) - diameter( 1 ) );
psi  = 0.5 * kr1 * diameter( 1 ) + psi2 - 0.5 * kr3 * diameter( 2 );
%  propagation functions
fac = k3 / k1 * sqrt( kr1 ./ kr3 ) .* exp( 1i * psi );
P.te = fac .* t12.te .* t23.te ./ ( 1 - r21.te .* r23.te .* exp( 2i * psi2 ) );
P.tm = fac .* t12.tm .* t23.tm ./ ( 1 - r21.tm .* r23.tm .* exp( 2i * psi2 ) );
%  propagate coefficients
a = P.te .* a;
b = P.tm .* b;
%  reshape coefficients
a = permute( reshape( a, siz( 2 ), siz( 1 ), [] ), [ 2, 1, 3 ] );
b = permute( reshape( b, siz( 2 ), siz( 1 ), [] ), [ 2, 1, 3 ] );


function [ r, t ] = fresnel( mat1, mat2, k0, kz )
%  FRESNEL - Fresnel coefficients.

%  wavevector in radial direction
kr1 = sqrt( mat1.k( k0 ) ^ 2 - kz .^ 2 );
kr2 = sqrt( mat2.k( k0 ) ^ 2 - kz .^ 2 );
%  permittivities and permeabilities
[ eps1, mu1 ] = deal( mat1.eps( k0 ), mat1.mu( k0 ) );
[ eps2, mu2 ] = deal( mat2.eps( k0 ), mat2.mu( k0 ) );
%  reflection coefficients, Hohenester Eq. (8.26)
r.te = (  mu2 * kr1 -  mu1 * kr2 ) ./ (  mu2 * kr1 +  mu1 * kr2 );
r.tm = ( eps2 * kr1 - eps1 * kr2 ) ./ ( eps2 * kr1 + eps1 * kr2 );
%  transmission coefficients
t.te = 2 *  mu2 * kr1 ./ (  mu2 * kr1 +  mu1 * kr2 ); 
t.tm = 2 * eps2 * kr1 ./ ( eps2 * kr1 + eps1 * kr2 ); 
