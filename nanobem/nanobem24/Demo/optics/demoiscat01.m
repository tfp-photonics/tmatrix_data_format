%  DEMOISCAT01 - iSCAT of nanoparticle above substrate.

%  material properties of layer structure and nanoparticles
mat1 = Material( 1.50 ^ 2, 1 );
mat2 = Material( 1.33 ^ 2, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
mat4 = Material( 1.59 ^ 2, 1 );
air  = Material( 1, 1 );
%  wavenumber of light in vacuum and glass
k0 = 2 * pi / 520;
k1 = mat1.k( k0 );

%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  nanosphere
diameter = 50;
p = trisphere( 144, diameter );
%  set switch value for different setups
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

%  planewave excitation
%    iSCAT if einc.dir( 3 ) points upwards
%    COBRI if einc.dir( 3 ) points downwards
pol = [ 1, 0, 0 ];
dir = [ 0, 0, 1 ];
%  excitation angle of 35°
einc = optics.roty( 35 ) * optics.decompose( k0, pol, dir );
%  solve BEM equation
bem = stratified.bemsolver( tau, layer, 'order', [] );
sol = bem \ einc( tau, 'layer', layer );

%  image lens, rotation matrix for optical axis in negative z-direction
NA = 1.3;
rot = optics.roty( 180 );
lens = optics.lensimage( mat1, air, k0, NA, 'rot', rot );
%  farfields and reflected incoming fields
far = farfields( sol, lens.dir );
refl = secondary( einc, layer, 'dir', 'down' );

% %  plot farfields
% plot( lens, far );

%  focus and image positions
x = 2000 * linspace( -1, 1, 201 );
%  image of fields for given focus position
focus = [ 0, 0, 0 ];
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
