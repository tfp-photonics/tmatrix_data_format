%  DEMOCAPILLARY02 - Same as democapillary01.m for various focus planes.

%  materials
mat1 = Material( 1.33 ^ 2, 1 );
mat2 = Material( 1.59 ^ 2, 1 );
mat3 = Material( 2.25, 1 );
mat4 = Material( 1, 1 );

%  polystrene nanosphere embedded in water
diameter = 1000;
%  wavenumber of light in vacuum
k0 = 2 * pi / 532;
%  Wiscombe cutoff for angular degree
x = 0.5 * diameter * mat1.k( k0 );
lmax = ceil( x + 2 + 4.05 * x ^ ( 1 / 3 ) );
%  T-matrix
mie = multipole.miesolver( mat2, mat1, diameter, lmax );
tmat = eval( mie, k0 );

%  Laguerre-Gauss beam with topological charge m=2, paraxial approximation
field = laguerregauss( mat1 );

field.foc = 50.8e6;
field.w0  = 1.8e6;
field.nLG = 0;
field.mLG = 2;
field.pow = 1.65;

%  shift for particle origin and function for incoming fields
shift = [ 5e3, 0, 0 ];
fun = @( pos, k0 ) paraxial( field, pos, k0, 'shift', shift );
%  solve Mie equations
q = qinc( tmat, fun, 'diameter', diameter );
sol = solve( tmat, q );

%  initialize FFT image lens 
rot = optics.roty( -90 );
lens = optics.lensimage2( mat4, mat4, k0, 0.4, 'rot', rot, 'n', 800 );
%  capillary 
mat = [ mat1, mat3, mat4 ];
diameter = [ 1.3, 1.8 ] * 1e6;
cap = optics.capillary( mat, diameter );
%  trace farfields through capillary
far = farfields( cap, sol, lens.dir, 'shift', shift );

%  focus length
z = ( 0 : 0.2 : 1 ) * 250e3;
%  loop over focus lenth
for it = 1 : 6
  %  image fields, shift focus 
  [ e, x ] = efield( lens, far, 'n', 2000, 'focus', [ z( it ), 0, 0 ] );  

  %  plot intensity
  subplot( 2, 3, it );
  imagesc( x * 1e-3, x * 1e-3, dot( e, e, 3 ) .' );

  set( gca, 'YDir', 'norm' );
  axis equal tight
  
  xlim( [ -150, 150 ] );
  ylim( [ -150, 150 ] );  

  xlabel( 'z (\mum)' );
  ylabel( 'y (\mum)' );
  title( [ 'z=', num2str( z( it ) * 1e-3 ), ' \mum' ] );
end
