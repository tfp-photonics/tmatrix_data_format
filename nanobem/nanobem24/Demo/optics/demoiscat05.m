%  DEMOISCAT05 - iSCAT images for different focus positions.

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
diameter = 55;
p = trisphere( 144, diameter );
mat = [ mat1, mat2, mat3 ];
p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  initialize BEM solver
bem = stratified.bemsolver( tau, layer, 'order', [] );
%  Jones matrix for quarter wavelength plate rotated by 45°
quarter = optics.quarterplate( 45 );
%  planewave excitation, apply quarter wave plate and rotation
einc = optics.decompose( k0, [ 1, 0, 0 ] * quarter .', [ 0, 0, 1 ] );
einc = optics.roty( 14 ) * einc;
%  solve BEM equation
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
x = linspace( - 2000, 2000, 101 );
z = linspace( - 4000, 4000, 51 );
%  allocate output
[ isca, iref ] = deal( zeros( numel( x ), numel( x ), 3, numel( z ) ) );

multiWaitbar( 'Focus', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over focus positions
for iz = 1 : numel( z )
  %  focus position
  focus = [ 0, 0, z( iz ) ];
  %  image of fields
  isca( :, :, :, iz ) = efield( lens, far,  x, x, 'focus', focus );
  iref( :, :, :, iz ) = efield( lens, refl, x, x, 'focus', focus );
    
  if ~mod( iz, 10 ),  multiWaitbar( 'Focus', iz / numel( z ) );  end
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%  flip x-axis because of rotation in lensimage
isca = flip( isca, 1 );
iref = flip( iref, 1 );
%  interference image
im = squeeze( dot( isca + iref, isca + iref, 3 ) - dot( iref, iref, 3 ) );

%%  plot iSCAT image together with slider
fig = uifigure;
g = uigridlayout( fig );
g.RowHeight = { '1x', 'fit' };
g.ColumnWidth = { '1x' };

ax = uiaxes( g );
h = imagesc( ax, 1e-3 * x, 1e-3 * x, im( :, :, 1 ) .' );

set( ax, 'YDir', 'norm' );
colormap( ax, gray( 200 ) );
 
title( ax, sprintf( 'focus = %5.2f mu', 1e-3 * z( 1 ) ) );
axis( ax, 'equal', 'tight' );
 
xlabel( ax, 'x (\mum)' );
ylabel( ax, 'y (\mum)' );

sld = uislider( g, 'Limits', [ 0, 50 ] );
sld.MajorTicks = [ 0, 10, 20, 30, 40, 50 ];
sld.ValueChangingFcn = @( src, event ) fun( src, event, ax, h, im, z );


function fun( ~, event, ax, h, im, z )
  %  FUN - Update function.
  iz = round( event.Value );
  set( h, 'CData', im( :, :, iz + 1 ) .' );
  title( ax, sprintf( 'focus = %5.2f mu', 1e-3 * z( iz + 1 ) ) );
  drawnow
end
