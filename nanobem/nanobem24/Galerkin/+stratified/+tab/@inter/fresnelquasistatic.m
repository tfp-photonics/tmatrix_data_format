function tq = fresnelquasistatic( obj, k0 )
%  FRESNELQUASISTATIC - Quasistatic Fresnel coefficients.

if abs( obj.i1 - obj.i2 ) == 1
  [ ~, tq ] = fresnel( obj.layer, k0, 1e10, obj.i2, obj.i1 );
else
  tq = struct( 'te', 0, 'tm', 0 );
end
