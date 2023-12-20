%  DEMOSTRATDIP01 - Decay rate for dipole and gold sphere above substrate.

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
%  dipole object
pos = bsxfun( @plus, mean( p.pos, 1 ), 0.7 * diameter * [ 1, 0, 0 ] );
pt = Point( mat, 2, pos );
dip = stratified.dipole( layer, pt );

%  light wavelength in vacuum
lambda = linspace( 400, 800, 20 );
k0 = 2 * pi ./ lambda;
%  allocate decay rate enhancement
tot = zeros( numel( k0 ), 3 );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavenumbers
for i = 1 : numel( k0 )
  %  solution of BEM equations
  sol = bem \ dip( tau, k0( i ) );
  %  total decay rate
  tot( i, : ) = decayrate( dip, sol );
 
  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );


%%  final plot
plot( lambda, tot, 'o-'  );  hold on

xlabel( 'Wavelength (nm)' );
ylabel( 'Total decay rate (P_0)' );
