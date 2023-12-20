function obj = fill2( obj, k0 )
%  FILL2 - Fill reflected Green function tables.
%
%  Usage for obj = stratified.tab.inter :
%    obj = fill2( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum

%  reflection coefficient for quasistatic approximation
obj.layer = eval( obj.layer, k0 );
tq = fresnelquasistatic( obj, k0 );
%  evaluate Green function table
[ r, z1, z2 ] = ndgrid( obj.rtab, obj.ztab1, obj.ztab2 );
y = eval( obj.somint, r, struct( 'z1', z1, 'z2', z2 ),  ...
  @( varargin ) fun( obj, tq, varargin{ : } ), k0, 'singular', 1 );

%  tables for interpolation
for name = fieldnames( y ) .'
  obj.ytab.( name{ 1 } ) =  ...
    griddedInterpolant( r, z1, z2, y.( name{ 1 } ), obj.method, 'nearest' );
end
%  save wavenumber
obj.k0 = k0;


function y = fun( obj, tq, data, kr, ~, mode )
  %  FUN - Integrand for Green function evaluation, Chew (6.15-16).

  %  unique radii
  [ r, ~, ind ] = unique( data.r );
  %  Bessel or Hankel functions
  switch mode
    case 'bessel'
      z0 = besselj( 0, kr * r );
      z1 = besselj( 1, kr * r );
    case 'hankel'
      z0 = besselh( 0, 1, kr * r );
      z1 = besselh( 1, 1, kr * r );
  end
  %  integrand w/o reflection coefficients
  z0 = reshape( z0( ind ), size( data.r ) );
  z1 = reshape( z1( ind ), size( data.r ) );
  %  secondary waves
  f = secondary2( obj, data, kr );
  
  %  propagation direction of reflected waves and propagation distance
  dir = sign( obj.i1 - obj.i2 );  
  Z = abs( data.z1 - data.z2 );
  %  quasistatic approximation
  q = stratified.zsqrt( - kr ^ 2 );
  f0q = 1i * exp( 1i * q * Z ) / ( 4 * pi ) .* z0;
  f1q = 1i * exp( 1i * q * Z ) / ( 4 * pi ) .* z1 * sign( real( kr ) );  
  
  %  Green function elements for SL
  y.te = f.te .* z0;
  y.tm = f.tm .* z0;
  %  first derivative along z, subtract quasistatic contribution
  y.tez1 = f.z1.te .* z0 - 1i * dir * tq.te * f0q;
  y.tmz1 = f.z1.tm .* z0 - 1i * dir * tq.tm * f0q;
  y.tez2 = f.z2.te .* z0 + 1i * dir * tq.te * f0q; 
  y.tmz2 = f.z2.tm .* z0 + 1i * dir * tq.tm * f0q;
  %  second derivative along z
  y.tezz = f.zz.te .* z0 - q * tq.te * f0q; 
  y.tmzz = f.zz.tm .* z0 - q * tq.tm * f0q;
  %  surface Green function, Chew (6.28)
  y.tes = kr ^ 2 * f.te .* z0 + q * tq.te * f0q;
  y.tms = kr ^ 2 * f.tm .* z0 + q * tq.tm * f0q;
  
  %  Green function elements for DL, radial derivative
  y.ter = - kr * f.te .* z1 - 1i * tq.te * f1q;
  y.tmr = - kr * f.tm .* z1 - 1i * tq.tm * f1q;
  %  mixed derivatives
  y.terz1 = - kr * f.z1.te .* z1 + dir * q * tq.te * f1q;
  y.terz2 = - kr * f.z2.te .* z1 - dir * q * tq.te * f1q;
  y.tmrz1 = - kr * f.z1.tm .* z1 + dir * q * tq.tm * f1q;
  y.tmrz2 = - kr * f.z2.tm .* z1 - dir * q * tq.tm * f1q;
  