function [ a, b ] = miecoefficients( obj, k0, varargin )
%  MIECOEFFICIENTS - Mie coefficients at sphere outside.
%
%  Usage for obj = feibelman.miesolver :
%    [ a, b ] = miecoefficients( obj, k0 )
%  Input
%    k0     :  wavevector of light in vacuum
%  Output
%    a,b    :  Mie coefficients for outside field

%  wavenumbers and permittivities
[ k1, eps1 ] = deal( obj.mat1.k( k0 ), obj.mat1.eps( k0 ) );
[ k2, eps2 ] = deal( obj.mat2.k( k0 ), obj.mat2.eps( k0 ) );
%  radius and Feibelman parameters
[ r, d1, d2 ] = deal( 0.5 * obj.diameter, obj.dperp, obj.dpar );
if isa( d1, 'function_handle'),  d1 = d1( k0 );  end
if isa( d2, 'function_handle'),  d2 = d2( k0 );  end

l = reshape( 1 : obj.lmax, [], 1 );
[ z1, z2 ] = deal( r * k1, r * k2 );
%  compute Riccati-Bessel functions
[ j1, zjp1 ] = riccatibessel( l, z1, 'j' );
[ j2, zjp2 ] = riccatibessel( l, z2, 'j' );
[ h2, zhp2 ] = riccatibessel( l, z2, 'h' );

%  correction terms for TM coefficients
r1 = ( eps1 - eps2 ) / r * ( j2 .* j1 .* l .* ( l + 1 ) * d1 + zjp2 .* zjp1 * d2 );
r2 = ( eps1 - eps2 ) / r * ( h2 .* j1 .* l .* ( l + 1 ) * d1 + zhp2 .* zjp1 * d2 );
%  Goncalves et al., Nat. Comm. 11, 366 (2020), Eq. (4a)
a = ( eps1 * j1 .* zjp2 - eps2 * j2 .* zjp1 + r1 ) ./   ...
    ( eps1 * j1 .* zhp2 - eps2 * h2 .* zjp1 + r2 );
  
%  correction terms for TE coefficients
r1 = ( z1 ^ 2 - z2 ^  2 ) / r * j2 .* j1 * d2;
r2 = ( z1 ^ 2 - z2 ^  2 ) / r * h2 .* j1 * d2;
%  Goncalves et al., Nat. Comm. 11, 366 (2020), Eq. (4b)
b = ( j1 .* zjp2 - j2 .* zjp1 + r1 ) ./ ( j1 .* zhp2 - h2 .* zjp1 + r2 );
