function rq = fresnelquasistatic( obj, k0 )
%  FRESNELQUASISTATIC - Quasistatic Fresnel coefficients.

rq = [ fresnel( obj.layer, k0, 1e10, obj.i1, obj.i1 - 1 ),  ...
       fresnel( obj.layer, k0, 1e10, obj.i1, obj.i1 + 1 ) ];
