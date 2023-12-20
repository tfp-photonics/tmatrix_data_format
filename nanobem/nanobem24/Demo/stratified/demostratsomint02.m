%  DEMOSOMINT02 - Sommerfeld integration, plot integrand w/o quasistatic.

%  material properties 
mat1 = Material( 2.25, 1 );
mat2 = Material( 1, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
%  wavenumber of light in vacuum
k0 = 2 * pi / 600;

layer = stratified.layerstructure( [ mat1, mat3, mat2 ], [ - 40, 0 ] );
%  Sommerfeld integrator
op = odeset( 'AbsTol', 1e-5, 'InitialStep', 1e-3 );
somint = stratified.isommerfeld( layer, 3, 'semi1', 5, 'op', op, 'cutoff', 0.1 );
    
r = 1e-3;
z = 0.1;

%  quasistatic Fresnel coefficient
rq = fresnel( layer, k0, 1e10, layer.n + 1, layer.n );
y = cumeval( somint, r, z,  ...
  @( varargin ) fun( layer, rq, varargin{ : } ), k0, 'singular', 1, 'nquad', 200 );

%  final plot
yt = horzcat( y.tmzz );
plot( abs( yt ), '.-' );  hold on

xlabel( 'Integration path' );
ylabel( 'Cumulative integral' );


function y = fun( layer, rq, data, kr, kz, mode )
  %  FUN - Integrand for scalar reflected Green function.
  switch mode
    case 'bessel'
      z0 = besselj( 0, kr * data.r );
    case 'hankel'
      z0 = besselh( 0, 1, kr * data.r );
  end
  %  quasistatic approximation
  q = stratified.zsqrt( - kr ^ 2 );
  f0 =  - 1i * kz * exp( 1i * q * data.z ) / ( 4 * pi ) * z0; 
  f0q = - 1i *  q * exp( 1i * q * data.z ) / ( 4 * pi ) * z0;  
  
  refl = rtcoeffs( layer, data.k0, kr, 'dir', 'down' );   
  y.tezz = refl.te * f0 - rq.te * f0q;
  y.tmzz = refl.tm * f0 - rq.tm * f0q;
end
