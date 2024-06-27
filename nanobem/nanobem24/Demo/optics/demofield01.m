%  DEMOFIELD01 - Field map for optically excited nanoparticle.

%  material properties of layer structure and nanoparticles
mat1 = Material( 1.50 ^ 2, 1 );
mat2 = Material( 1.33 ^ 2, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
mat4 = Material( 1.59 ^ 2, 1 );
%  wavenumber of light in vacuum and glass
k0 = 2 * pi / 520;
k1 = mat1.k( k0 );

%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  nanosphere
diameter = 50;
p = trisphere( 144, diameter );

switch 1
  case 1
    %  gold nanosphere 5 nm above substrate
    mat = [ mat1, mat2, mat3 ];
    p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
    %  boundary elements with linear shape functions
    tau = BoundaryEdge( mat, p, [ 3, 2 ] );
  case 2
    %  coated gold nanosphere 5 nm above substrate
    mat = [ mat1, mat2, mat3, mat4 ];
    p2 = trisphere( 144, diameter + 10 );
    p. verts( :, 3 ) = p. verts( :, 3 ) + 0.5 * diameter + 10;
    p2.verts( :, 3 ) = p2.verts( :, 3 ) + 0.5 * diameter + 10;
    %  flip faces for proper simulation of coated sphere
    p.faces = fliplr( p.faces );
    tau = BoundaryEdge( mat, p2, [ 4, 2 ], p, [ 4, 3 ] );    
  case 3
    %  PS nanosphere 5 nm above substrate
    mat = [ mat1, mat2, mat4 ];
    p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
    %  boundary elements with linear shape functions
    tau = BoundaryEdge( mat, p, [ 3, 2 ] );  
  case 4
    %  coupled gold nanospheres 5 nm above substrate
    mat = [ mat1, mat2, mat3, mat3 ];
    p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
    [ p1, p2 ] = deal( p );
    %  gap distance
    gap = 10;
    p1.verts( :, 1 ) = p1.verts( :, 1 ) - 0.5 * diameter - 0.5 * gap;
    p2.verts( :, 1 ) = p2.verts( :, 1 ) + 0.5 * diameter + 0.5 * gap;
    %  boundary elements with linear shape functions
    tau = BoundaryEdge( mat, p1, [ 3, 2 ], p2, [ 4, 2 ] );     
end

%  initialize BEM solver
bem = stratified.bemsolver( tau, layer, 'order', [] );
%  planewave excitation
%    iSCAT if dir( 3 ) points upwards
%    COBRI if dir( 3 ) points downwards
t = 20 / 180 * pi;
pol = [ - cos( t ), 0, sin( t ) ];
dir = [   sin( t ), 0, cos( t ) ];
einc = optics.decompose( k0, pol, dir );
%  solve BEM equation
sol = bem \ einc( tau, 'layer', layer );


%  grid in xz-plane
x = linspace( -100, 100, 51 );
z = linspace( -100, 150, 71 );
[ xx, zz ] = ndgrid( x, z );
%  place points in dielectric environment
pos = [ xx( : ), 0 * xx( : ), zz( : ) ];
pt = stratified.Point( layer, tau, pos );
%  comoute electromagnetic fields
[ e1, h1 ] = fields( sol, pt );
[ e2, h2 ] = fields( einc, pt, 'layer', layer );

%  plot electric field
figure
ee = reshape( real( e1( :, 1 ) + e2( :, 1 ) ), size( xx ) );

imagesc( x, z, ee .' );  hold on
plot( [ min( x ), max( x ) ], [ 0, 0 ], 'w-' );

colormap redblue( 200 );
set( gca, 'YDir', 'norm' );
axis equal tight

xlabel( 'x (nm)' );
ylabel( 'z (nm)' );
