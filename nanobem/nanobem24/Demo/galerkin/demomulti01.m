%  DEMOMULTI01 - T-matrix for dielectric nanosphere and single wavelength.

%  material properties
mat1 = Material( 1, 1 );
mat2 = Material( 4, 1 );
%  material vector
mat = [ mat1, mat2 ];

%  nanosphere
%    when needed use different number of vertices, e.g., 144, 400, 900
diameter = 160;
n = 144;
p = trisphere( n, diameter );
%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 2, 1 ] );     
 %  wavenumber of light in vacuum
k0 = 2 * pi / 1000;
%  BEM solver and T-matrix solver
lmax = 4;
bem = galerkin.bemsolver( tau, 'order', [] );
tsolver = multipole.tsolver( mat, 1, lmax );
%  T-matrix
sol = bem \ tsolver( tau, k0 );
t1 = eval( tsolver, sol );
%par = particle( p.verts, p.faces );
%savegmsh( par, "mesh.msh");
%  additional information for H5 file
info = multipole.h5info(tau);
info.name = "Sphere test";
info.keywords = "";
info.description = "Single sphere and single wavelength";
info.matname = [ "Air, Vacuum", "Custom" ];
info.shape = "sphere";
info.params.radius = diameter/2; 
info.unit = "nm";
info.method_parameters.number_points = n;
%  save T-matrix
h5save( t1, 'ctmatrix_sphere.h5', info );


%  comparison with Mie theory
mie = multipole.miesolver( mat2, mat1, diameter, lmax );
t2 = eval( mie, k0 );

%  full T-matrices
r1 = full( t1 );
r2 = full( t2 );
fun = @abs;
%  final plot
plot( fun( diag( r1 ) ), 'o' ); hold on
plot( fun( diag( r2 ) ), '+' );

legend( 'T-matrix', 'Mie' );

xlabel( '#Mode' );
ylabel( 'T_{nn}' );

set( gca, 'YScale', 'log' );
