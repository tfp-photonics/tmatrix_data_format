function obj = fill2( obj, k0 )
%  FILL2 - Fill reflected Green function tables.
%
%  Usage for obj = stratified.tab.intra1 :
%    obj = fill2( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum

%  reflection coefficient for quasistatic approximation
obj.layer = eval( obj.layer, k0 );
rq = fresnelquasistatic( obj, k0 );
%  evaluate Green function table
[ r, z ] = ndgrid( obj.rtab, obj.ztab );
y = eval( obj.somint, r, z,  ...
  @( varargin ) fun( obj, rq, varargin{ : } ), k0, 'singular', 1 );

%  tables for interpolation
for name = fieldnames( y ) .'
  obj.ytab.( name{ 1 } ) =  ...
    griddedInterpolant( r, z, y.( name{ 1 } ), obj.method, 'nearest' );
end
%  save wavenumber
obj.k0 = k0;



function y = fun( obj, rq, data, kr, kz, mode )
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
  f0 = 1i * exp( 1i * kz * data.z ) / ( 4 * pi * kz ) .* z0;
  f1 = 1i * exp( 1i * kz * data.z ) / ( 4 * pi * kz ) .* z1;
  %  quasistatic approximation
  q = stratified.zsqrt( - kr ^ 2 );
  f0q = 1i * exp( 1i * q * data.z ) / ( 4 * pi ) .* z0;
  f1q = 1i * exp( 1i * q * data.z ) / ( 4 * pi ) .* z1 * sign( real( kr ) );  
  %  generalized reflection coefficients
  r = secondary( obj.layer, data.k0, kr, obj.i1, obj.i1 );
  switch obj.i1
    case 1
      dir = - 1;
    otherwise
      dir = + 1;
  end
  
  %  Green function elements for SL
  y.te = r.te * f0;
  y.tm = r.tm * f0;
  %  first derivative along z, subtract quasistatic contribution
  y.tez = dir * ( 1i * kz * r.te * f0 - 1i * rq.te * f0q );
  y.tmz = dir * ( 1i * kz * r.tm * f0 - 1i * rq.tm * f0q );
  %  second derivative along z
  y.tezz = - kz ^ 2 * r.te * f0 + q * rq.te * f0q;
  y.tmzz = - kz ^ 2 * r.tm * f0 + q * rq.tm * f0q;
  %  surface Green function, Chew (6.28)
  y.tes = kr ^ 2 * r.te * f0 + q * rq.te * f0q;
  y.tms = kr ^ 2 * r.tm * f0 + q * rq.tm * f0q;
   
  %  Green function elements for DL, radial derivative
  y.ter = - kr * r.te * f1 - 1i * rq.te * f1q;
  y.tmr = - kr * r.tm * f1 - 1i * rq.tm * f1q;
  %  mixed derivatives
  y.terz = dir * ( - 1i * kz * kr * r.te * f1 + q * rq.te * f1q );
  y.tmrz = dir * ( - 1i * kz * kr * r.tm * f1 + q * rq.tm * f1q );
