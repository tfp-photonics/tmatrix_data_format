function a = amplitude( obj, k0 )
%  AMPLITUDE - Field amplitude for given input power and light frequency.
%
%  Usage for obj = laguerregauss :
%    a = amplitude( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum
%  Output
%    a    :  field amplitude

%  focal parameter from calculated waist
k = obj.mat.k( k0 );
f = 1 / ( k * waist( obj, k0 ) );
%  unnormalized energy flux
flux = 1 / k ^ 2 * integral( @( t ) obj.fun( sin( t ) /  ...
  ( sqrt( 2 ) * f ) ) .^ 2 .* ( 1 + cos( t ) .^ 2 ) .* sin( t ), 0, 0.5 * pi );
%  compute amplitude
Z = 376.730 * obj.mat.Z( k0 );
a = sqrt( 2 * Z * obj.pow / ( pi * flux ) );
