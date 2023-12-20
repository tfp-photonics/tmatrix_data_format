function obj = fill2( obj, k0 )
%  FILL2 - Fill reflected Green function tables.
%
%  Usage for obj = stratified.tab.intra2 :
%    obj = fill2( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum

%  reflection coefficient for quasistatic approximation
obj.layer = eval( obj.layer, k0 );
rq = fresnelquasistatic( obj, k0 );
%  evaluate Green function tables
if all( obj.ztab1 == obj.ztab2 )
  %  evaluate Green function table
  [ r, z1 ] = ndgrid( obj.rtab, obj.ztab1 );
  z2 = z1;
  y = eval( obj.somint, r, z1,  ...
    @( varargin ) fun( obj, rq, [ 1, 2 ], varargin{ : } ), k0, 'singular', 1 );
  %  disentangle Y1 and Y2
  for name = fieldnames( y ) .'
    switch name{ 1 }( end )
      case '1'
        y1.( name{ 1 } ) = y.( name{ 1 } );
      case '2'
        y2.( name{ 1 } ) = y.( name{ 1 } );
    end
  end
else
  %  evaluate Green function table
  [ r, z1 ] = ndgrid( obj.rtab, obj.ztab1 );
  [ ~, z2 ] = ndgrid( obj.rtab, obj.ztab2 );
  y1 = eval( obj.somint, r, z1,  ...
      @( varargin ) fun( obj, rq, 1, varargin{ : } ), k0, 'singular', 1 );
  y2 = eval( obj.somint, r, z2,  ...
      @( varargin ) fun( obj, rq, 2, varargin{ : } ), k0, 'singular', 1 );
end

%  tables for interpolation
for name = fieldnames( y1 ) .'
  obj.ytab1.( name{ 1 }( 1 : end - 1 ) ) =  ...
    griddedInterpolant( r, z1, y1.( name{ 1 } ), obj.method, 'nearest' );
end
for name = fieldnames( y2 ) .'
  obj.ytab2.( name{ 1 }( 1 : end - 1 ) ) =  ...
    griddedInterpolant( r, z2, y2.( name{ 1 } ), obj.method, 'nearest' );  
end
%  save wavenumber
obj.k0 = k0;



function y = fun( obj, rq, it, data, kr, kz, mode )
  %  FUN - Integrand for Green function evaluation.

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
  %  reshape Bessel or Hankel functions
  z0 = reshape( z0( ind ), size( data.r ) );
  z1 = reshape( z1( ind ), size( data.r ) );  
  %  layer thickness
  layer = obj.layer;
  d = layer.d( obj.i1 - 1 );
  %  up- and downgoing waves
  f1 = 1i * exp( 1i * kz * data.z ) / ( 4 * pi * kz );
  f2 = 1i * exp( 1i * kz * ( 2 * d - data.z ) ) / ( 4 * pi * kz );
  
  %  quasistatic approximation
  q = stratified.zsqrt( - kr ^ 2 );
  f1q = 1i * exp( 1i * q * data.z ) / ( 4 * pi );
  f2q = 1i * exp( 1i * q * ( 2 * d - data.z ) ) / ( 4 * pi );
  sig = sign( real( kr ) );  
   
  %  secondary waves
  layer = obj.layer;
  wave = secondary( layer, data.k0, kr, obj.i1, obj.i1 );
  [ te, tm ] = deal( wave.te, wave.tm );
  
  %  z1 - z2
  if any( it == 1 )
    %  Green function elements for SL
    y.te1 = ( te( 1, 1 ) * f1 + te( 2, 2 ) * f2 ) .* z0;
    y.tm1 = ( tm( 1, 1 ) * f1 + tm( 2, 2 ) * f2 ) .* z0;
    %  first derivative along z
    y.tez1 = 1i * kz * ( te( 1, 1 ) * f1 - te( 2, 2 ) * f2 ) .* z0;
    y.tmz1 = 1i * kz * ( tm( 1, 1 ) * f1 - tm( 2, 2 ) * f2 ) .* z0;
    %  second derivative along z
    y.tezz1 = - kz ^ 2 * ( te( 1, 1 ) * f1 + te( 2, 2 ) * f2 ) .* z0;
    y.tmzz1 = - kz ^ 2 * ( tm( 1, 1 ) * f1 + tm( 2, 2 ) * f2 ) .* z0;
    %  surface Green function, Chew (6.28)
    y.tes1 = kr ^ 2 * ( te( 1, 1 ) * f1 + te( 2, 2 ) * f2 ) .* z0;
    y.tms1 = kr ^ 2 * ( tm( 1, 1 ) * f1 + tm( 2, 2 ) * f2 ) .* z0;
   
    %  Green function elements for DL, radial derivative
    y.ter1 = - kr * ( te( 1, 1 ) * f1 + te( 2, 2 ) * f2 ) .* z1;
    y.tmr1 = - kr * ( tm( 1, 1 ) * f1 + tm( 2, 2 ) * f2 ) .* z1;
    %  mixed derivatives
    y.terz1 = - 1i * kz * kr * ( te( 1, 1 ) * f1 - te( 2, 2 ) * f2 ) .* z1;
    y.tmrz1 = - 1i * kz * kr * ( tm( 1, 1 ) * f1 - tm( 2, 2 ) * f2 ) .* z1;
  end
  
  %  z1 + z2
  if any( it == 2 )
    %  Green function elements for SL
    y.te2 = ( te( 1, 2 ) * f1 + te( 2, 1 ) * f2 ) .* z0;
    y.tm2 = ( tm( 1, 2 ) * f1 + tm( 2, 1 ) * f2 ) .* z0;  
    %  first derivative along z, subtract quasistatic contribution
    y.tez2 = 1i * kz * ( te( 1, 2 ) * f1 - te( 2, 1 ) * f2 ) .* z0 -  ...
                1i * ( rq( 1 ).te * f1q - rq( 2 ).te * f2q ) .* z0;
    y.tmz2 = 1i * kz * ( tm( 1, 2 ) * f1 - tm( 2, 1 ) * f2 ) .* z0 -  ...
                1i * ( rq( 1 ).tm * f1q - rq( 2 ).tm * f2q ) .* z0; 
    %  second derivative along z
    y.tezz2 = - kz ^ 2 * ( te( 1, 2 ) * f1 + te( 2, 1 ) * f2 ) .* z0 +  ...
                   q * ( rq( 1 ).te * f1q + rq( 2 ).te * f2q ) .* z0;
    y.tmzz2 = - kz ^ 2 * ( tm( 1, 2 ) * f1 + tm( 2, 1 ) * f2 ) .* z0 +  ...
                   q * ( rq( 1 ).tm * f1q + rq( 2 ).tm * f2q ) .* z0;
    %  surface Green function, Chew (6.28)
    y.tes2 = kr ^ 2 * ( te( 1, 2 ) * f1 + te( 2, 1 ) * f2 ) .* z0 +  ...
                q * ( rq( 1 ).te * f1q + rq( 2 ).te * f2q ) .* z0;   
    y.tms2 = kr ^ 2 * ( tm( 1, 2 ) * f1 + tm( 2, 1 ) * f2 ) .* z0 +  ...
                q * ( rq( 1 ).tm * f1q + rq( 2 ).tm * f2q ) .* z0;                
  
    %  Green function elements for DL, radial derivative
    y.ter2 = - kr * ( te( 1, 2 ) * f1 + te( 2, 1 ) * f2 ) .* z1 -  ...
       1i * sig * ( rq( 1 ).te * f1q + rq( 2 ).te * f2q ) .* z1;
    y.tmr2 = - kr * ( tm( 1, 2 ) * f1 + tm( 2, 1 ) * f2 ) .* z1 -  ...
       1i * sig * ( rq( 1 ).tm * f1q + rq( 2 ).tm * f2q ) .* z1;
    %  mixed derivatives
    y.terz2 = - 1i * kz * kr * ( te( 1, 2 ) * f1 - te( 2, 1 ) * f2 ) .* z1 +  ...
                   q * sig * ( rq( 1 ).te * f1q - rq( 2 ).te * f2q ) .* z1;
    y.tmrz2 = - 1i * kz * kr * ( tm( 1, 2 ) * f1 - tm( 2, 1 ) * f2 ) .* z1 +  ...
                   q * sig * ( rq( 1 ).tm * f1q - rq( 2 ).tm * f2q ) .* z1;
  end    
    