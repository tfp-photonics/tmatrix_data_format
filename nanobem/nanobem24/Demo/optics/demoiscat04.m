%  DEMOISCAT04 - iSCAT of nanosphere above substrate w/ wave-quarter plate.

%  material properties of layer structure and nanoparticles
mat1 = Material( 1.50 ^ 2, 1 );
mat2 = Material( 1.33 ^ 2, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
air  = Material( 1, 1 );
%  wavenumber of light in vacuum and glass
k0 = 2 * pi / 520;
k1 = mat1.k( k0 );

%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  gold nanosphere 5 nm above substrate
diameter = 50;
p = trisphere( 144, diameter );
p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  Jones matrix for quarter-wavelength plate rotated by 45°
quarter = optics.quarterplate( 45 );
%  planewave excitation
einc = optics.decompose( k0, [ 1, 0, 0 ] * quarter .', [ 0, 0, 1 ] );
einc = optics.roty( 20 ) * einc;    %  rotate incoming field
%  solve BEM equation
bem = stratified.bemsolver( tau, layer, 'order', [] );
sol = bem \ einc( tau, 'layer', layer );

%  image lens and rotation matrix, optical axis in negative z-direction
NA = 1.3;
rot = optics.roty( 180 );
lens = optics.lensimage( mat1, air, k0, NA, 'rot', rot, 'mcut', 5 );
lens.backfocal = @( x ) x * quarter .';
%  farfields and reflected incoming fields
far = farfields( sol, lens.dir );
refl = secondary( einc, layer, 'dir', 'down' );

%  focus and image positions
x = 2000 * linspace( -1, 1, 201 );
z = 0;
%  focus position
focus = [ 0, 0, z ];
%  image of fields
isca = efield( lens, far,  x, x, 'focus', focus );
iref = efield( lens, refl, x, x, 'focus', focus );
%  flip x-axis because of rotation in lensimage
isca = flip( isca, 1 );
iref = flip( iref, 1 );
%  interference image
im = dot( isca + iref, isca + iref, 3 );


%  final figure
figure
imagesc( 1e-3 * x, 1e-3 * x, im .' );

set( gca, 'YDir', 'norm' );
axis equal tight
colormap gray( 200 )

xlabel( 'x (\mum)' );
ylabel( 'y (\mum)' );
