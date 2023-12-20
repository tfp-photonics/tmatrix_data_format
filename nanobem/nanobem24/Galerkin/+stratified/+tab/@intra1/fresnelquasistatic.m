function rq = fresnelquasistatic( obj, k0 )
%  FRESNELQUASISTATIC - Quasistatic Fresnel coefficients.

switch obj.i1
  case 1
    rq = fresnel( obj.layer, k0, 1e10, 1 );
  otherwise
    n = obj.layer.n;
    rq = fresnel( obj.layer, k0, 1e10, n + 1, n ); 
end
