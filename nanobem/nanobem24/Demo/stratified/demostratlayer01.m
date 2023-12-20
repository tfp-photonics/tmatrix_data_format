%  DEMOSTRATLAYER01 - Secondary fields in layer structure.

%  material properties of layer structure
mat1 = Material( 2.25, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );
mat3 = Material( 1, 1 );
%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat, [ 0, 40 ] );

%  wavenumber of light in vacuum and parallel wavenumber
k0 = 2 * pi / 600;
kpar = 0.4 * k0;
%  field and source positions
z1 = linspace( - 500, 800, 501 );
z2 = - 100;
%  primary and secondary fields
f = fields( layer, k0, kpar, z1, z2, 'primary', 1 );

%  plot TM fields
figure
plot( z1, real( f.te ), '.-', z1, real( f.tm ), '.-' );  hold on

legend( 'TE', 'TM' );
xlim( [ min( z1 ), max( z1 ) ] );

xlabel( 'z_1 (nm)' );
ylabel( 'Field' );


