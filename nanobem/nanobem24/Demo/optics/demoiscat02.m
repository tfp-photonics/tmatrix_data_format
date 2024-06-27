%  DEMOISCAT02 - iSCAT of nanosphere, change focus of imaging lens.

%  material properties of layer structure andnanoparticles
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

%  planewave excitation
%    iSCAT if dir( 3 ) points upwards
%    COBRI if dir( 3 ) points downwards
einc = optics.roty( 0 ) * optics.decompose( k0, [ 1, 0, 0 ], [ 0, 0, 1 ] );
%  solve BEM equation
bem = stratified.bemsolver( tau, layer, 'order', [] );
sol = bem \ einc( tau, 'layer', layer );

%  image lens, rotation matrix for optical axis in negative z-direction
NA = 1.3;
rot = optics.roty( 180 );
lens = optics.lensimage( mat1, air, k0, NA, 'rot', rot, 'mcut', 5 );
%  farfields and reflected incoming fields
far = farfields( sol, lens.dir );
refl = secondary( einc, layer, 'dir', 'down' );

%  focus and image positions
x = 2000 * linspace( -1, 1, 201 );
y = 0;
z = linspace( - 4000, 4000, 201 );
%  allocate output
[ isca, iref ] = deal( zeros( numel( x ), numel( y ), 3, numel( z ) ) );

multiWaitbar( 'Focus', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over focus positions
for iz = 1 : numel( z )
  %  focus position
  focus =  [ 0, 0, z( iz ) ];
  %  image of fields
  isca( :, :, :, iz ) = efield( lens, far,  x, y, 'focus', focus );
  iref( :, :, :, iz ) = efield( lens, refl, x, y, 'focus', focus );
    
  if ~mod( iz, 10 ),  multiWaitbar( 'Focus', iz / numel( z ) );  end
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%  flip x-axis because of rotation in lensimage
isca = flip( isca, 1 );
iref = flip( iref, 1 );
%  interference image and interference phase 
im1 = dot( isca + iref, isca + iref, 3 ) - dot( iref, iref, 3 );
im2 = angle( dot( isca, iref, 3 ) );

%  interfernce image
figure
i1 = ceil( 0.5 * numel( y ) );
imagesc( 1e-3 * x, z * 1e-3, squeeze( im1( :, i1, : ) ) .' );

set( gca, 'YDir', 'norm' );
axis equal tight
colormap gray( 200 )
colorbar( 'northoutside' );

xlabel( 'x (\mum)' );
ylabel( 'z (\mum)' );
set( gca, 'FontSize', 11 );

%  phase angle image
figure
i1 = ceil( 0.5 * numel( y ) );
imagesc( 1e-3 * x, z * 1e-3, squeeze( im2( :, i1, : ) ) .' );

set( gca, 'YDir', 'norm' );
axis equal tight
mycolormap( 'cen:3' );
h = colorbar( 'northoutside' );
set( h, 'Ticks', [ - 0.99 * pi, 0, 0.99 * pi ] );
set( h, 'TickLabels', { '-\pi', '0', '\pi' } );
set( h, 'FontSize', 11 );

xlabel( 'x (\mum)' );
set( gca, 'YTickLabels', [], 'FontSize', 11 );
% ylabel( 'z (\mum)' );
