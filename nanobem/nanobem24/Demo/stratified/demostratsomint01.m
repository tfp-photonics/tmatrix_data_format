%  DEMOSTRATSOMINT01 - Sommerfeld integration, plot integrand.

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

% y = eval( somint, r, z,  ...
%   @( varargin ) fun( layer, varargin{ : } ), k0, 'singular', 1 )

y = cumeval( somint, r, z,  ...
  @( varargin ) fun( layer, varargin{ : } ), k0, 'singular', 1, 'nquad', 200 );

%  final plot
yt = horzcat( y.tmzz );
semilogy( abs( yt ), '.-' );  hold on

xlabel( 'Integration path' );
ylabel( 'Cumulative integral' );


function y = fun( layer, data, kr, kz, mode )
  %  FUN - Integrand for scalar reflected Green function.
  switch mode
    case 'bessel'
      z0 = besselj( 0, kr * data.r );
    case 'hankel'
      z0 = besselh( 0, 1, kr * data.r ); 
  end
  refl = rtcoeffs( layer, data.k0, kr, 'dir', 'down' );   
  y.tezz = - 1i * kz * refl.te * z0 * exp( 1i * kz * data.z ) / ( 4 * pi );
  y.tmzz = - 1i * kz * refl.tm * z0 * exp( 1i * kz * data.z ) / ( 4 * pi );
end
