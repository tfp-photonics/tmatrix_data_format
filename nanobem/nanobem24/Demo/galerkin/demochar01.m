%  DEMOCIM01 - Characteristic modes for optically excited ellipsoid.

%  nanoellipsoid
diameter = 50;
p = trisphere( 144, diameter );
p.verts( :, 1 ) = p.verts( :, 1 ) * 2;

mat1 = Material( 1, 1 );
mat2 = Material( epsdrude( 'Au' ), 1 );
%  material properties
mat = [ mat1, mat2 ];

%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 2, 1 ] );
%  BEM and characteristic mode solver
bem = galerkin.bemsolver( tau );
modes = charModes( bem );

%  planewave excitation
exc = galerkin.planewave( [ 1, 0, 0 ], [ 0, 0, 1 ] );

%  light wavelength in vacuum
lambda = linspace( 400, 700, 20 );
k0 = 2 * pi ./ lambda;
%  allocate optical cross sections, number of characteristic modes
[ csca1, csca2, num ] = deal( zeros( 1, numel( k0 ) ) );
%  tolerance for cutoff of characteristic modes
tol = 1e-2;

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavenumbers
for i = 1 : numel( k0 )
  %  evaluate characteristic modes and solve BEM equations
  modes = eval( modes, k0( i ) );
  sol1 = solve( modes, exc( tau, k0( i ) ) );
  %  expansion coefficients for characteristic modes
  a = expand( modes, sol1 );
  %  project on largest modes
  ind = abs( a ) >= tol * max( abs( a ) );
  sol2 = project( modes, sol1, ind );
  num( i ) = nnz( ind );
  
  %  optical cross sections
  csca1( i ) = scattering( exc, sol1 );
  csca2( i ) = scattering( exc, sol2 );
  
  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%%  final plot
figure
plot( lambda, csca1, 'o-' );  hold on
plot( lambda, csca2, '+-' );

legend( 'BEM', 'charModes' );

xlabel( 'Wavelength (nm)' );
ylabel( 'Cross section (nm^2)' );
