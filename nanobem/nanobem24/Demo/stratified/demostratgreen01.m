%  DEMOSTRATGREEN01 - Tabulated Green function for sphere above substrate.

%  material properties of layer structure
mat1 = Material( 2.25, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );
mat3 = Material( 1, 1 );
%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  nanosphere
p = trisphere( 144, 50 );
p.verts( :, 3 ) = p.verts( :, 3 ) - min( p.verts( :, 3 ) ) + 5;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  slice positions and set up computational grid
pos = vertcat( tau.pos );
r = slice( layer, pos, pos );
r = grid( layer, r, 'nr', 20, 'nz', 20 );
%  set up tabulated refletced Green function object
green = stratified.tab.green( r );

%  wavenumber of light in vacuum 
k0 = 2 * pi / 600;
%  fill Green function and plotting
green = fill( green, k0 );
plot( green );
