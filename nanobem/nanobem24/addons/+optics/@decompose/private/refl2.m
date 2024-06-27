function [ r, t ] = refl2( layer, k0, dir, einc )
%  REFL1 - Reflected and transmitted fields for downgoing waves.
%
%  Usage :
%    [ r, t ] = refl1( layer, k0, dir, einc )
%  Input
%    layer    :  layer structure
%    k0       :  wavenumber of light in vacuum
%    dir      :  wave propagation direction
%    einc     :  wave amplitude
%  Output
%    r,t      :  reflected and transmitted waves

%  normalization function
norm = @( x ) bsxfun( @rdivide, x, sqrt( dot( x, x, 2 ) ) );
%  triad of propagation direction and TE, TM vectors
ez = repmat( [ 0, 0, 1 ], size( dir, 1 ), 1 );
te = norm( cross( dir, ez, 2 ) );
if any( isnan( te( : ) ) )
  [ te( isnan( te( :, 1 ) ), : ) ] = deal( [ 1, 0, 0 ] );
end
tm = norm( cross( dir, te, 2 ) ); 

%  materials and interface positions in lowest and uppermost medium
[ mat1, z1 ] = deal( layer.mat(   1 ), layer.z(   1 ) );
[ mat2, z2 ] = deal( layer.mat( end ), layer.z( end ) );
%  wavenumbers and impedances
[ k1, Z1 ] = deal( mat1.k( k0 ), mat1.Z( k0 ) );
[ k2, Z2 ] = deal( mat2.k( k0 ), mat2.Z( k0 ) );
%  parallel and z-components of wavevector
kpar = k2 * hypot( dir( :, 1 ), dir( :, 2 ) );
k1z = stratified.zsqrt( k1 ^ 2 - kpar .^ 2 );
k2z = stratified.zsqrt( k2 ^ 2 - kpar .^ 2 );
%  generalized transmission and reflection coefficients
[ r, t ] = rtcoeffs( layer, k0, kpar, 'dir', 'down' ); 

%  propagation directions for reflected waves
dir1 = [ dir( :, 1 : 2 ), - dir( :, 3 ) ];
%  reflected fields
tm1 = - cross( dir1, cross( dir, tm ) );
e1 = bsxfun( @times, te,  r.te .* dot( te, einc, 2 ) ) +  ...
     bsxfun( @times, tm1, r.tm .* dot( tm, einc, 2 ) );
e1 = bsxfun( @times, e1, exp( - 2i * k2z * z2 ) );   

%  propagation direction in lowest medium
dir2 = [ k2 * dir( :, 1 : 2 ), - k1z ] / k1;
i2 = ~imag( dir2( :, 3 ) );
%  transmitted fields
tm2 = - cross( dir2, cross( dir, tm ) );
e2 = bsxfun( @times, te,  t.te .* dot( te, einc, 2 ) ) +  ...
     bsxfun( @times, tm2, t.tm .* dot( tm, einc, 2 ) ) * Z1 / Z2;
e2 = bsxfun( @times, e2, exp( - 1i * ( k2z * z2 - k1z * z1 ) ) );

%  set output, keep only propagating waves
r = struct( 'efield', e1, 'dir', dir1 );
t = struct( 'efield', e2( i2, : ), 'dir', dir2( i2, : ) );
