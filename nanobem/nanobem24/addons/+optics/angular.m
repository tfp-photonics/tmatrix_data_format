function far = angular( mat, k0, dir, e, x, y, varargin )
%  ANGULAR - Farfields from electric fields in given plane.
%
%  Usage :
%    far = optics.angular( mat, k0, dir, e, x, y, PropertyPairs )
%  Input
%    mat      :  material properties
%    k0       :  wavenumber of light in vacuum      
%    dir      :  farfield propagation directions
%    e        :  electric field in plane
%    x,y      :  real-space grid
%  PropertyName
%    z        :  vertical position of plane
%    nmax     :  slice number of directions into bunches of size NMAX
%  Output
%    far      :  farfields in direction of DIR

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'z', 0 );
addParameter( p, 'nmax', 1000 );
%  parse input
parse( p, varargin{ : } );

%  wavenumber in medium and allocate output
k = mat.k( k0 );
far = zeros( size( dir ) );
%  increments for grid
hx = x( 2 ) - x( 1 );
hy = y( 2 ) - y( 1 );
%  real space grid
[ x, y ] = ndgrid( x, y );
[ x, y, e ] = deal( x( : ), y( : ), reshape( e, [], 3 ) );
%  initialize iterator
it = 0;

while it < size( dir, 1 )
  %  index to directions
  ind = it + ( 1 : p.Results.nmax );
  ind = ind( ind <= size( dir, 1 ) );
  %  Fourier transform of electric field
  fac = k * dir( ind, 1 ) * x .' + k * dir( ind, 2 ) * y .';
  far( ind, : ) = exp( - 1i * fac ) * e;
  %  update iterator
  it = it + numel( ind );
end

%  wavevector in z-direction
kz = k * dir( :, 3 );
%  add prefactors and convert Fourier transform to farfield
far = far * hx * hy / ( 4 * pi ^ 2 );
far = - 2i * pi *  ...
  bsxfun( @times, far, kz .* exp( - 1i * kz * p.Results.z ) );
  