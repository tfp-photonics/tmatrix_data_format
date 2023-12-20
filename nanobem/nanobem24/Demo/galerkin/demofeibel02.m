%  DEMOFEIBEL02 - Electric field of optically excited coupled nanospheres.

%  diameter of sphere and gap distance
diameter = 20;
gap = 5;
%  spheres
[ p1, p2 ] = deal( trisphere( 144, diameter ) );
p1.verts( :, 3 ) =   p1.verts( :, 3 ) - 0.5 * ( diameter + gap );
p2.verts( :, 3 ) = - p2.verts( :, 3 ) + 0.5 * ( diameter + gap );
%  flip faces 
p2.faces = p2.faces( :, [ 1, 3, 2 ] );

%  dielectric functions
mat1 = Material( 1.33 ^ 2, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );
mat3 = Material( 4, 1 );
%  material properties
mat = [ mat1, mat2, mat3 ];

%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p1, [ 2, 1 ], p2, [ 3, 1 ] );
%  index to boundary elements of Au nanosphere
inout = vertcat( tau.inout );
ind = inout( :, 1 ) == 2;
%  constant Feibelman parameters for Au only
[ d1, d2 ] = deal( - 0.4, 0 );
param = feibelman.param( tau( ind ), @( ~ ) deal( d1, d2 ) );
%  initialize BEM solver
bem = feibelman.bemsolver( tau, param, 'relcutoff', 2, 'order', [] );
%  planewave excitation
exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

% %  solve BEM equation
k0 = 2 * pi / 600;
sol = bem \ exc( tau, k0 );

%  points
n = 41;
xx = 0.8 * diameter * linspace( - 1, 1,     n );
zz = 0.8 * diameter * linspace( - 2, 2, 2 * n );
[ x, z ] = ndgrid( xx, zz );
pt = Point( tau, [ x( : ), 0 * x( : ), z( : ) ] );
%  evaluate fields
[ e, h ] = fields( sol, pt, 'relcutoff', 2, 'waitbar', 1 );
[ ei, hi ] = fields( exc, pt, k0 );

%%
%  final plot
ee = real( squeeze( dot( e + ei, e + ei, 2 ) ) );
%  plot intensity
figure
imagesc( xx, zz, reshape( ee, size( x ) ) .' );  hold on

t = linspace( 0, 2 * pi, 101 );
z0 = 0.5 * diameter + 0.5 * gap;
plot( 0.5 * diameter * cos( t ), 0.5 * diameter * sin( t ) - z0, 'w-' );
plot( 0.5 * diameter * cos( t ), 0.5 * diameter * sin( t ) + z0, 'w-' );

set( gca, 'YDir', 'norm' );
axis equal tight

xlabel( 'x (nm)' );
ylabel( 'z (nm)' );
