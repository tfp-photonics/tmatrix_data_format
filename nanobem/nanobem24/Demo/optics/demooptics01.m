%  DEMOOPTICS01 - Focused laser beam.

%  material properties 
mat1 = Material( 1.33 ^ 2, 1 );
%  wavenumber of light in vacuum
k0 = 2 * pi / 520;

%  focus lens
NA = 0.8;
lens = optics.lensfocus( mat1, k0, NA );
%  incoming fields, before crossing the Gaussian reference sphere
e = normpdf( lens.rho, 0, 0.8 );
e = e( : ) * [ 1, 0, 0 ];

%  points where fields are computed
x = 2000 * linspace( -1, 1, 101 );
z = 4000 * linspace( -1, 1, 201 );
[ xx, zz ] = ndgrid( x, z );
pts = Point( mat1, 1, [ xx( : ), 0 * xx( : ), zz( : ) ] );

%%  focal fields
figure
foc = eval( lens, e );
e1 = fields( foc, pts );
%  plot focal fields
subplot( 1, 3, 1 );
imagesc( x, z, reshape( dot( e1, e1, 2 ), size( xx ) ) .' );

set( gca, 'YDir', 'norm' );
axis equal tight

%%  user-defined focus position
foc = eval( lens, e, 'focus', [ 0, 0, 1000 ] );
e1 = fields( foc, pts );
%  plot focal fields
subplot( 1, 3, 2 );
imagesc( x, z, reshape( dot( e1, e1, 2 ), size( xx ) ) .' );  hold on
plot( 0, 1000, 'm+' );

set( gca, 'YDir', 'norm' );
axis equal tight

%%  rotate focal fields
lens.rot = optics.roty( 30 );
foc = eval( lens, e );
e1 = fields( foc, pts );
%  plot focal fields
subplot( 1, 3, 3 );
imagesc( x, z, reshape( dot( e1, e1, 2 ), size( xx ) ) .' );

set( gca, 'YDir', 'norm' );
axis equal tight
