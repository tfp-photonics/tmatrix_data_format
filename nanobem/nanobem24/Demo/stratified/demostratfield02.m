%  DEMOSTRATFIELD02 - Electric field for optically excited nanodisk.

%  material properties of layer structure
mat1 = Material( 2.25, 1 );
mat2 = Material( 1, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  polygon for disk
poly = polygon( 25, 'size', [ 80, 80 ] ); 
%  edge profile for nanodisk
edge = edgeprofile( 30, 10, 'mode', '11' );  
%  extrude polygon to nanoparticle
p = tripolygon( poly, edge );
p.verts( :, 3 ) = p.verts( :, 3 ) - min( p.verts( :, 3 ) ) + 1e-2;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  initialize BEM solver
bem = stratified.bemsolver( tau, layer, 'waitbar', 1 );
%  planewave excitation
exc = stratified.planewave( layer, [ 1, 0, 0 ], [ 0, 0, -1 ] );
%  light wavelength in vacuum
lambda = 620;
k0 = 2 * pi / lambda;

%  solution of BEM equation
sol = bem \ exc( tau, k0 );

%  positions where field is computed
x = linspace( - 70, 70,  51 );
z = linspace( - 50, 100, 51 );
[ xx, zz ] = ndgrid( x, z );
pos = [ xx( : ), 0 * xx( : ), zz( : ) ];
pt = stratified.Point( layer, tau, pos );
%  comoute electromagnetic fields
[ e1, h1 ] = fields( sol, pt );
[ e2, h2 ] = fields( exc, pt, k0 );

%%
%  plot electric field intensity
ee = reshape( dot( e1 + e2, e1 + e2, 2 ), size( xx ) );

imagesc( x, z, log10( ee ) .' );  hold on
plot( [ - 70, 70 ], [ 0, 0 ], 'w-' );

set( gca, 'YDir', 'norm' );
axis equal tight

xlabel( 'x (nm)' );
ylabel( 'z (nm)' );