%  DEMOISCAT06 - iSCAT of gold nanosphere above substrate.
%    See Mahmoodabadi et al., Opt. Express 28, 25969 (2020), Fig. 3(b).

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

%  nanosphere
diameter = 50;
p = trisphere( 144, diameter );
%  gold nanosphere 3000 nm above substrate
mat = [ mat1, mat2, mat3 ];
p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 3000;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  planewave excitation
pol = [ 1, 0, 0 ];
dir = [ 0, 0, 1 ];
einc = optics.decompose( k0, pol, dir );
%  solve BEM equation
bem = stratified.bemsolver( tau, layer, 'order', [] );
sol = bem \ einc( tau, 'layer', layer );

%  image lens, rotation matrix for optical axis in negative z-direction
NA = 1.4;
rot = optics.roty( 180 );
lens = optics.lensimage( mat1, air, k0, NA, 'rot', rot );
%  farfields and reflected incoming fields
far = farfields( sol, lens.dir );
refl = secondary( einc, layer, 'dir', 'down' );

%  focus and image positions
x = 2000 * linspace( -1, 1, 201 );
y = 0;
z = linspace( 0, 10000, 201 );
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
%  interference image
im = dot( isca + iref, isca + iref, 3 );

%  interfernce image
figure
i1 = ceil( 0.5 * numel( y ) );
imagesc( 1e-3 * x, z * 1e-3, squeeze( im( :, i1, : ) ) .' );  hold on
plot( 1e-3 * [ min( x ), max( x ) ], [ 3, 3 ], 'k--' );

set( gca, 'YDir', 'norm' );
axis equal tight

colormap gray( 200 )
colormap( flipud( colormap ) );
colorbar( 'eastoutside' );

xlabel( 'x (\mum)' );
ylabel( 'z (\mum)' );
set( gca, 'FontSize', 11 );
