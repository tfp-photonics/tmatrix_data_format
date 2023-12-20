%  DEMOFEIBEL01 - Optical cross section for sphere w/ Feibelman parameters.

%  nanosphere
diameter = 20;
p = trisphere( 144, diameter );

%  dielectric functions for water and gold
mat1 = Material( 1.33 ^ 2, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );
%  material properties
mat = [ mat1, mat2 ];

%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 2, 1 ] );

%  constant Feibelman parameters
[ d1, d2 ] = deal( - 0.4, 0 );
param = feibelman.param( tau, @( ~ ) deal( d1, d2 ) );
%  BEM solver
rules = quadboundary.rules( 'quad3', triquad( 1 ) );
bem = feibelman.bemsolver( tau, param, 'rules', rules, 'waitbar', 1 );

%  light wavelength in vacuum
lambda = linspace( 400, 800, 20 );
k0 = 2 * pi ./ lambda;
%  planewave excitation
exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

%  allocate optical cross section
cext = deal( zeros( numel( k0 ), size( exc.pol, 1 ) ) );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavenumbers
for i = 1 : numel( k0 )
  %  solve BEM equations
  sol = bem \ exc( tau, k0( i ) );
  %  exctinction cross section
  cext( i, : ) = extinction( exc, sol );
  
  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%%  final plot
plot( lambda, cext, 'o-'  );  hold on

xlabel( 'Wavelength (nm)' );
ylabel( 'Cross section (nm^2)' );

%  set up Mie solver with Feibelman parameters
mie = feibelman.miesolver( mat2, mat1, diameter, 'lmax', 20 );
mie.dperp = d1;
mie.dpar  = d2;
%  set up Mie solver w/o Feibelman parameters
mie2 = miesolver( mat2, mat1, diameter );

plot( lambda, extinction( mie, k0 ), '+-' );
set( gca, 'ColorOrderIndex', 2 );
plot( lambda, extinction( mie2, k0 ), 'o--' );
legend( 'BEM', 'Mie', 'Mie^{(0)}' );
