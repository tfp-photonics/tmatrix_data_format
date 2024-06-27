%  DEMOFOCUS01 - Focusing lens and field imaging.

%  material properties 
mat1 = Material( 1.5 ^ 2, 1 );
mat2 = Material( 1, 1 );
%  wavenumber of light in vacuum
k0 = 2 * pi / 520;

%  focus lens
NA = 0.4;
lens = optics.lensfocus( mat1, k0, NA );
%  incoming fields
e = normpdf( lens.rho, 0, 0.5 ) .* exp( 1i * lens.phi );
e = e( : ) * [ 1, 0, 0 ];
%  planewave decomposition of focal fields
foc = eval( lens, e );

%  points in focal plane
x = 2000 * linspace( -1, 1, 101 );
[ xx, yy ] = ndgrid( x );
pts = Point( [ mat1, mat2 ], 1, [ xx( : ), yy( : ), 0 * xx( : ) ] );
%  evaluate focal fields
[ e1, h1 ] = fields( foc, pts );

%  lens for imaging
NA = 0.6;
lens2 = optics.lensimage( mat1, mat2, k0, NA );
%  image farfields from angular spectrum representation and fields from
%  planewave decomposition
far = optics.angular( mat1, k0, lens2.dir, e1, x, x, 0 );
i1 = efield( lens2, far, x, x );
i2 = efield( lens2, foc, x, x );

%  Poynting vector in z-direction, compare z=0 and image plane
P1 = 0.5 * real( cross( e1, conj( h1 ), 2 ) );
P1 = sum( P1( :, 3 ) );     %  ideally P1 should equal P2
P2 = 0.5 * mat2.Z( k0 ) * sum( reshape( dot( i1, i1, 3 ), [], 1 ) );

%  final plot
figure

for data = struct( 'it', { 1, 2 }, 'e', { i1, i2 } )
  %  plot intensity
  subplot( 1, 2, data.it ); 
  imagesc( 1e-3 * x, 1e-3 * x, dot( data.e, data.e, 3 ) .' );
  colorbar

  set( gca, 'YDir', 'norm' );
  axis equal tight

  xlabel( 'x (\mum)' );
  ylabel( 'y (\mum)' );
end

