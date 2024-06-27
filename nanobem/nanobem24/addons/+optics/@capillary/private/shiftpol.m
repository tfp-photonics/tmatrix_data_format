function [ a2, b2, M2 ] = shiftpol( a1, b1, M1, t, x )
%  SHIFTPOL - Shift origin for polar scattering coefficients.
%
%  Usage :
%    [ a2, b2, M2 ] = shiftpol( a1, b1, M1, t, x )
%  Input
%    a1,b1  :  polar scattering coefficients
%    M1     :  angular orders
%    t      :  polar angles
%    x      :  k * shift
%  Output
%    a2,b2  :  polar scattering coefficients for shifted origin
%    M2     :  angular orders for shifted origin

%  angular orders for shift operation, Abramowitz (9.3.1)
xmax = max( hypot( x( :, 1 ), x( :, 2 ) ) ) * max( sin( t ) );
m1 = max( M1 );
m2 = max( m1, ceil( 1.3591 * xmax ) );
M2 = reshape( - m2 : m2, [], 1 ); 
%  allocate output
nz = numel( a1 ) / ( numel( M1 ) * numel( t ) );
[ a2, b2 ] = deal( zeros( numel( M2 ), numel( t ), nz ) );

for iz = 1 : nz
  %  shift in z-direction
  a = bsxfun( @times, a1( :, :, iz ), exp( - 1i * x( iz, 3 ) * cos( t .' ) ) );
  b = bsxfun( @times, b1( :, :, iz ), exp( - 1i * x( iz, 3 ) * cos( t .' ) ) );
  %  convert transverse shift parameter to polar coordinates
  [ u, r ] = cart2pol( x( iz, 1 ) * sin( t ), x( iz, 2 ) * sin( t ) );
  
  %  Fourier coefficients for shift operator
  shift = besselj2( M2, r ) .* exp( - 1i * M2 * u .' );
  n = numel( M2 ) + numel( M1 );
  %  convolution
  a = ifft( fft( shift, n, 1 ) .* fft( a, n, 1 ) );
  b = ifft( fft( shift, n, 1 ) .* fft( b, n, 1 ) );
  %  extract output
  a2( :, :, iz ) = a( m1 + 1 : end - m1 - 1, : );
  b2( :, :, iz ) = b( m1 + 1 : end - m1 - 1, : );
end
 