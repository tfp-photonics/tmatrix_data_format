%  DEMOOPTICS02 - Image of dipoles, optical axis +z.

%  material properties 
mat = Material( 1.33 ^ 2, 1 );
air = Material( 1, 1 );
%  wavenumber of light in vacuum
k0 = 2 * pi / 520;

%  dipoles
dip1 = galerkin.dipole( Point( mat, 1, [ -500, -100, 0 ] ) );
dip2 = galerkin.dipole( Point( mat, 1, [  500,  100, 1000 ] ) );
%  imaging lens, image dipoles from above
NA = 0.8;
lens = optics.lensimage( mat, air, k0, NA );
dir = lens.dir;

%  far-fields of dipole
pt = Point( mat, 1, dir );
far1 = squeeze( farfields( dip1, pt, k0 ) );
far2 = squeeze( farfields( dip2, pt, k0 ) );
far = far1( :, :, 1 ) + far2( :, :, 2 );    %  dipoles oriented along x,y

%  image positions, focus planes
n = 101;
x = 1500 * linspace( -1, 1, n );
zf = -500 : 500 : 1000;

%%  image of dipoles
figure

%  scan focus plane through planes where dipoles are located
for iz = 1 : numel( zf )
  %  image of dipoles
  e = efield( lens, far, x, x, 'focus', [ 0, 0, zf( iz ) ] );
  %  plot image of dipoles
  subplot( 2, 2, iz );
  imagesc( x, x, dot( e, e, 3 ) .' );
  
  set( gca, 'YDir', 'norm' );
  axis equal tight
  title( [ 'zf=', num2str( zf( iz ) ), ' nm' ] );
end
