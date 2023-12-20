%  DEMOSTRATLAYER02 - Plane wave propagating through layer structure.

%  material properties of layer structure
mat1 = Material( 2.25, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );
mat3 = Material( 1, 1 );
%  set up layer structure
mat = [ mat1, mat2, mat3 ];
layer = stratified.layerstructure( mat, [ 0, 40 ] );

%  wavenumber of light in vacuum 
k0 = 2 * pi / 600;
%  fields
z = linspace( - 500, 800, 501 );
pos = z( : ) * [ 0, 0, 1 ];
%  plane waves inside of layerstructure
pol = [ 1, 0, 0 ];
dir = [ 0, 0, -1 ];
[ e, h ] = planewave( layer, pol, dir, k0, pos );

%  plot TE and TM fields
figure
plot( z, real( e( :, 1 ) ), '.-', z, real( h( :, 2 ) ), '.-' );  

legend( 'TE', 'TM' );
xlim( [ min( z ), max( z ) ] );

xlabel( 'z (nm)' );
ylabel( 'Field' );