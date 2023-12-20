function yout = cumeval2( obj, r, z, fun, k0, nquad )
%  CUMEVAL2 - Cumulative Sommerfeld integral w/ singular value subtraction.
%
%  Usage for obj = stratified.isommerfeld :
%    yout = cumeval2( obj, r, z, fun, k0, nquad )
%  Input
%    r,z        :  radii and z-values for evaluation 
%    fun        :  Sommerfeld integrand
%    k0         :  wavenumber of light in vacuum
%    nquad      :  number of quadrature points
%  Output
%    yout       :  cumulative Sommerfeld integral

%  wavenumbers 
k = arrayfun( @( mat ) mat.k( k0 ), obj.layer.mat, 'uniform', 1 );
[ k, kmax ] = deal( k( obj.i1 ), obj.semi1 * max( real( k ) ) );
%  integrand for kr = 0
data = select( obj, r, z, k0, k, 'semi' );
yout = fun( data, 0, k, 'bessel' );
y0 = reshape( compress( obj, yout ), numel( r ), [] );

%  integration to real or imaginary infinity
[ data1, ind1 ] = select( obj, r, z, k0, k, 'real' );
[ data2, ind2 ] = select( obj, r, z, k0, k, 'imag' );
%  integration functions
fun1 = @( x, ~ ) odefun2( obj, data,  x, kmax, fun, y0, 'semi' );
fun2 = @( x, ~ ) odefun2( obj, data1, x, kmax, fun, y0( ind1, : ), 'bessel' );
fun3 = @( x, ~ ) odefun2( obj, data2, x, kmax, fun, y0( ind2, : ), 'hankel' );

%  perform integration along semi-ellipse
x = linspace( 1e-10, pi, nquad );
[ ~, yt ] = ode45( fun1, x, log( obj.cutoff * k0 ) * y0( : ), obj.op );
%  reshape output
yt = reshape( yt, nquad, numel( r ), [] );
y0 = reshape( yt( end, :, : ), numel( r ), [] );

%  perform integration to positive infinity
x = fliplr( logspace( - 10, 0, nquad ) );
if nnz( ind1 )
  [ ~, yt1 ] = ode45( fun2, x, y0(  ind1, : ), obj.op );
  yt( nquad + 1 : 2 * nquad, ind1, : ) = reshape( yt1, nquad, nnz( ind1 ), [] );
end
%  perform integration to imaginary infinity
if nnz( ind2 )
  [ ~, yt2 ] = ode45( fun3, x, y0( ind2, : ), obj.op );
  yt( nquad + 1 : 2 * nquad, ind2, : ) = reshape( yt2, nquad, nnz( ind2 ), [] );
end

%  set output
yt = mat2cell( yt, ones( 1, 2 * nquad ), size( yt, 2 ), size( yt, 3 ) );
yout = cellfun( @( y ) uncompress( obj, y, yout ), yt, 'uniform', 1 );

