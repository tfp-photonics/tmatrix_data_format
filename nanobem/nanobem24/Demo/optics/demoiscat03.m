%  DEMOISCAT03 - iSCAT of gold nanosphere, comparison with dipole model.

%  material properties of layer structure and gold nanosphere
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
einc = optics.roty( 20 ) * optics.decompose( k0, [ 1, 0, 0 ], [ 0, 0, 1 ] );
%  solve BEM equation
bem = stratified.bemsolver( tau, layer, 'order', [] );
sol = bem \ einc( tau, 'layer', layer );

%  electric field at sphere center
pt = Point( mat, 2, mean( p.pos, 1 ) );
edip = fields( einc, pt, 'layer', layer );
%  polarizability of sphere, Hohenester Eq. (9.9)
[ epsb, epsz ] = deal( mat2.eps( k0 ), mat3.eps( k0 ) );
alpha = 0.5 * pi * diameter ^ 3 * epsb * ( epsz - epsb ) / ( epsz + 2 * epsb );
%  dipole moment of sphere
P = alpha * edip;

% %  dipole moment from surface charge
% P = sum( p.area .* surfc( sol ) .* p.pos, 1 );

%  dipole object
dip = stratified.dipole( layer, pt );

%  image lens, rotation matrix for optical axis in negative z-direction
NA = 1.3;
rot = optics.roty( 180 );
lens = optics.lensimage( mat1, Material( 1, 1 ), k0, NA, 'rot', rot );
%  farfields and reflected incoming fields
far1 = farfields( sol, lens.dir );
far2 = farfields( dip, lens.dir, k0 );
far2 = reshape( reshape( far2, [], 3 ) * P( : ), [], 3 );
refl = secondary( einc, layer, 'dir', 'down' );

%  focus and image positions
x = 2000 * linspace( -1, 1, 201 );
z = 0;
%  focus position
focus = [ 0, 0, z ];
%  image of fields, flip x-axis because of rotation in lensimage
isca1 = flip( efield( lens, far1, x, x, 'focus', focus ), 1 );
isca2 = flip( efield( lens, far2, x, x, 'focus', focus ), 1 );
iref  = flip( efield( lens, refl, x, x, 'focus', focus ), 1 );
%  interference image
im1 = dot( isca1 + iref, isca1 + iref, 3 );
im2 = dot( isca2 + iref, isca2 + iref, 3 );

%  final figure
figure
name = [ "BEM", "Dipole" ];
for data = struct( 'it', { 1, 2 }, 'im', { im1, im2 } )
  %  plot image
  subplot( 1, 2, data.it );
  imagesc( 1e-3 * x, 1e-3 * x, data.im .' );

  set( gca, 'YDir', 'norm' );
  axis equal tight
  caxis( [ min( [ im1( : ); im2( : ) ] ),  ...
           max( [ im1( : ); im2( : ) ] ) ] );
  colorbar

  xlabel( 'x (\mum)' );
  ylabel( 'y (\mum)' );
  title( name( data.it ) );
end
%  load colormap
colormap gray( 200 )

% %  plot cut through image
% figure
% plot( x, im1( :, 101 ), x, im2( :, 101 ) );
% 
% legend( 'BEM', 'Dipole' );
% 
% xlabel( 'x (\mum' );
% ylabel( 'Intensity' );
