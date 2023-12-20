function yout = eval1( obj, r, z, fun, k0 )
%  EVAL1 - Evaluate Sommerfeld integral w/o singular value subtraction.
%
%  Usage for obj = stratified.isommerfeld :
%    yout = eval1( obj, r, z, fun, k0 )
%  Input
%    r,z        :  radii and z-values for evaluation 
%    fun        :  Sommerfeld integrand
%    k0         :  wavenumber of light in vacuum
%  Output
%    yout       :  evaluated Sommerfeld integral

%  wavenumbers 
k = arrayfun( @( mat ) mat.k( k0 ), obj.layer.mat, 'uniform', 1 );
[ k, kmax ] = deal( k( obj.i1 ), obj.semi1 * max( real( k ) ) );
%  integrand for kr = 0
data = select( obj, r, z, k0, k, 'semi' );
yout = fun( data, 0, k, 'bessel' );

%  integration to real or imaginary infinity
[ data1, ind1 ] = select( obj, r, z, k0, k, 'real' );
[ data2, ind2 ] = select( obj, r, z, k0, k, 'imag' );
%  integration functions
fun1 = @( x, ~ ) odefun1( obj, data,  x, kmax, fun, 'semi' );
fun2 = @( x, ~ ) odefun1( obj, data1, x, kmax, fun, 'bessel' );
fun3 = @( x, ~ ) odefun1( obj, data2, x, kmax, fun, 'hankel' );

%  perform integration along semi-ellipse
y0 = 0 * compress( obj, yout );
[ ~, y ] = ode45( fun1, [ 1e-10, 1e-3, pi ], y0( : ), obj.op );
y = reshape( y( end, : ), numel( r ), [] );
%  perform integration to positive infinity
if nnz( ind1 )
  [ ~, y1 ] = ode45( fun2, [ 1, 1e-3, 1e-10 ], y( ind1, : ), obj.op );
  y( ind1, : ) = reshape( y1( end, : ), nnz( ind1 ), [] );
end
%  perform integration to imaginary infinity
if nnz( ind2 )
  [ ~, y2 ] = ode45( fun3, [ 1, 1e-3, 1e-10 ], y( ind2, : ), obj.op );
  y( ind2, : ) = reshape( y2( end, : ), nnz( ind2 ), [] );
end

%  set output
yout = uncompress( obj, y, yout );

