function f = secondary2( obj, data, kr )
%  SECONDARY2 - Secondary waves for reflected Green function.
%
%  Usage for obj = stratified.tab.inter :
%    wave = secondary2( obj, data, kr )
%  Input
%    data   :  auxiliary information from Sommerfeld integrator
%    kr     :  wavenumber in radial direction
%  Output
%    f      :  secondary TE and TM waves

%  coefficients for secondary waves
layer = obj.layer;
wave = secondary( layer, data.k0, kr, obj.i1, obj.i2 ); 
%  wavenumbers 
kz1 = wave.kz( obj.i1 );
kz2 = wave.kz( obj.i2 );
%  allocate output
[ f.te, f.tm, f.z1.te, f.z1.tm, f.z2.te,  ...
              f.z2.tm, f.zz.te, f.zz.tm ] = deal( 0 * data.r );
%  wave direction
dir = [ 1, -1 ];

%  loop over TE and TM polarizations
%    see also stratified.layerstructure/efield
for name = [ "te", "tm" ]
  a = wave.( name );
  %  loop over wave directions
  for d1 = 1 : size( a, 1 )
  for d2 = 1 : size( a, 2 )
    dir1 = dir( d1 );  if obj.i1 == 1,            dir1 = -1;  end
    dir2 = dir( d2 );  if obj.i2 == layer.n + 1,  dir2 = -1;  end
    %  distance to interfaces
    Z1 = fun( obj, data.z1, obj.i1,   dir1 );
    Z2 = fun( obj, data.z2, obj.i2, - dir2 );
    %  add up secondary fields
    fac = 1i * a( d1, d2 ) *  ...
      exp( 1i * ( kz1 * Z1 + kz2 * Z2 ) ) / ( 4 * pi * kz2 );
    f.( name ) = f.( name ) + fac; 
    %  z-derivatives
    f.z1.( name ) = f.z1.( name ) + 1i * kz1 * dir1 * fac;
    f.z2.( name ) = f.z2.( name ) - 1i * kz2 * dir2 * fac;
    f.zz.( name ) = f.zz.( name ) + kz1 * dir1 * kz2 * dir2 * fac;
  end
  end
end


function z = fun( obj, z, it, dir )
%  FUN - Distance from interfaces.
switch dir
  case 1
    z = z - obj.layer.z( it - 1 );
  otherwise
    z = obj.layer.z( it ) - z; 
end
