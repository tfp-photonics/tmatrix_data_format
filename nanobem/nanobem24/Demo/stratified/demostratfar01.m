%  DEMOSTRATFAR01 - Electric farfields for optically excited nanosphere.

%  material properties of layer structure
mat1 = Material( 2.25, 1 );
mat2 = Material( 1, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

%  nanosphere
diameter = 50;
p = trisphere( 144, diameter );
p.verts( :, 3 ) = p.verts( :, 3 ) + 0.5 * diameter + 5;
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 3, 2 ] );

%  initialize BEM solver
bem = stratified.bemsolver( tau, layer, 'waitbar', 1 );
%  planewave excitation
exc = stratified.planewave( layer, [ 1, 0, 0 ], [ 0, 0, -1 ] );
%  light wavelength in vacuum
lambda = 550;
k0 = 2 * pi / lambda;

%  solution of BEM equation
sol = bem \ exc( tau, k0 );

%  directions where electromagnetic farfields are computed
t = reshape( linspace( 0, 2 * pi, 201 ), [], 1 );
dir = [ sin( t ), 0 * t, cos( t ) ];
%  electromagnetic farfields
[ e, h ] = farfields( sol, dir );
%  Poynting vector in farfield direction
P = dot( dir, 0.5 * real( cross( e, conj( h ), 2 ) ), 2 );

%  polar plot
polarplot( t + 0.5 * pi, P );
