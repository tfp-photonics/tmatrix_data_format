%  DEMOCORR01 - Fit best focus plane to experimental data.

%  load experimental data
%    imExp    -  experimental iSCAT image
%    xExp     -  pixel positions in nm
C = load( 'dataExp' );

%  material properties of layer structure and nanoparticles
mat1 = Material( 1.50 ^ 2, 1 );
mat2 = Material( 1.00 ^ 2, 1 );
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
%  Jones matrix for quarter-wavelength plate rotated by 45°
quarter = optics.quarterplate( 45 );
%  planewave excitation, apply quarter-wave plate and rotation
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
x = C.xExp;
y = C.xExp;
z = linspace( - 4000, 4000, 101 );
%  allocate output
[ isca, iref ] = deal( zeros( numel( x ), numel( y ), 3, numel( z ) ) );

multiWaitbar( 'Focus', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over focus positions
for iz = 1 : numel( z )
  %  focus position
  focus = [ 0, 0, z( iz ) ];
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
%  interference image
im = squeeze( dot( isca + iref, isca + iref, 3 ) - dot( iref, iref, 3 ) );

%  correlation map
corr = zeros( 1, size( im, 3 ) );
%  compare simulated and experimental images
for it = 1 : size( corr, 2 )
  corr( it ) = corr2( im( :, :, it ), C.imExp );
end
%  maximal correlation
[ ~, i1 ] = max( corr );

%  final figure
figure
subplot( 1, 2, 1 );
imagesc( 1e-3 * x, 1e-3 * y, im( :, :, i1 ) .' );  

set( gca, 'YDir', 'norm' );
axis equal tight

xlabel( 'x (\mum)' );
ylabel( 'y (\mum)' );
title( 'Simulation' );

subplot( 1, 2, 2 );
imagesc( 1e-3 * x, 1e-3 * y, C.imExp .' );  

set( gca, 'YDir', 'norm' );
axis equal tight

xlabel( 'x (\mum)' );
ylabel( 'y (\mum)' );
title( 'Experiment' );

colormap gray( 200 )
