function y = odefun2( obj, data, x, kmax, fun, y0, ipath )
%  ODEFUN2 - ODE function with singular value subtraction.
%
%  Usage for obj = stratified.isommerfeld :
%    y = odefun1( obj, data, x, kmax, fun, ipath )
%  Input
%    data   :  auxiliary data for integration
%    x      :  integration variable
%    kmax   :  maximal wavenumber of layer structure
%    fun    :  user-defined integrand
%    y0     :  integrand for kr = 0
%    ipath  :  integration path 'semi', 'bessel', or 'hankel'
%  Output
%    y      : integrand

%  cutoff parameter for singularity subtraction
a = obj.cutoff * data.k0;

switch ipath
  case 'semi'
    %  integration along semi-ellipse
    kr = kmax * ( 1 - cos( x ) - 1i * obj.semi2 * sin( x ) );
    w = kmax * ( sin( x ) - 1i * obj.semi2 * cos( x ) );
    %  evaluate integrand
    kz = stratified.zsqrt( data.k ^ 2 - kr ^ 2 );
    y = w / kr * compress( obj, fun( data, kr, kz, 'bessel' ) );
    %  singularity subtraction
    y = y - w / kr * a ^ 2 / ( a ^ 2 + kr ^ 2 ) * y0( : );
  case 'bessel'
    %  integration to real infinity
    kr = 2 * kmax / x;
    w = - 2 * kmax / x ^ 2;
    %  evaluate integrand
    kz = stratified.zsqrt( data.k ^ 2 - kr ^ 2 );
    y = w / kr * compress( obj, fun( data, kr, kz, 'bessel' ) );
    %  singularity subtraction
    y = y - w / kr * a ^ 2 / ( a ^ 2 + kr ^ 2 ) * y0( : );    
  case 'hankel'
    %  integration to imaginary infinity
    kr1 = 2 * kmax * (   1 - 1i + 1i / x );
    kr2 = 2 * kmax * ( - 1 - 1i + 1i / x );
    w = - 2i * kmax / x ^ 2;
    %  evaluate integrand, Hohenester Eq. (B.10)
    kz1 = stratified.zsqrt( data.k ^ 2 - kr1 .^ 2 );
    kz2 = stratified.zsqrt( data.k ^ 2 - kr2 .^ 2 ); 
    y = 0.5 * w / kr1 * compress( obj, fun( data, kr1, kz1, 'hankel' ) ) -  ...
        0.5 * w / kr2 * compress( obj, fun( data, kr2, kz2, 'hankel' ) ) ;
    %  singularity subtraction
    y = y - w / kr1 * a ^ 2 / ( a ^ 2 + kr1 ^ 2 ) * y0( : ) +  ...
            w / kr2 * a ^ 2 / ( a ^ 2 + kr2 ^ 2 ) * y0( : );
end
