%  DEMOSTRATSPEC01 - Spectrum for gold sphere above glass substrate.

%  material properties of layer structure
mat1 = Material( 2.25, 1 );
mat2 = Material( 1, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  nanosphere
diameter = 50;
p = trisphere( 144, diameter );
p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  initialize BEM solver
bem = stratified.bemsolver( tau, layer, 'waitbar', 1 );
%  planewave excitation
exc = stratified.planewave( layer, [ 1, 0, 0 ], [ 0, 0, -1 ] );
%  light wavelength in vacuum
lambda = linspace( 400, 800, 20 );
k0 = 2 * pi ./ lambda;
%  allocate optical cross sections
[ csca, cabs ] = deal( zeros( size( k0 ) ) );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavenumbers
for i = 1 : numel( k0 )
  %  solution of BEM equations
  sol = bem \ exc( tau, k0( i ) );
  %  optical cross sections
  csca( i ) = scattering( exc, sol );
  cabs( i ) = absorption( exc, sol );
 
  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );


%%  final plot
plot( lambda, csca, 'o-'  );  hold on

xlabel( 'Wavelength (nm)' );
ylabel( 'Scattering cross section (nm^2)' );

%  comparison with Mie theory
mie = miesolver( mat3, mat2, diameter );
csca0 = scattering( mie, k0 );
cabs0 = absorption( mie, k0 );

plot( lambda, csca0, 's-' ); 
legend( 'BEM w/ layer', 'Mie w/o layer' );
