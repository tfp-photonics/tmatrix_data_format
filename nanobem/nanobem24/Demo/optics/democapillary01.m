%  DEMOCAPILLARY01 - Side image of nanoparticle through capillary.
%  See Hohenester et al., Nanophotonics 13, 457 (2024).

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

%  image fields, shift focus 
[ e, x ] = efield( lens, far, 'n', 2000, 'focus', [ 150e3, 0, 0 ] );  

%  plot intensity
figure
imagesc( x * 1e-3, x * 1e-3, dot( e, e, 3 ) .' );

set( gca, 'YDir', 'norm' );
axis equal tight

xlim( [ -150, 150 ] );
ylim( [ -150, 150 ] );

xlabel( 'z (\mum)' );
ylabel( 'y (\mum)' );
