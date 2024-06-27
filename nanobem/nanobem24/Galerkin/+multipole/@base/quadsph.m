function quad = quadsph( obj, key )
%  QUADSPH - Quadrature points and weights for spherical harmonics.
%
%  Usage for obj = multipole.base :
%    quad = quadsphere( obj, key )
%  Input
%    key      :  'matrix' (default) or 'vector' for output format
%  Output
%    quad.u   :  azimuthal angles
%    quad.t   :  polar angles
%    quad.w   :  quadrature weights

%  maximal angular order
tab = obj.tab;
lmax = max( tab.l ); 
%  quadrature points for azimuthal direction
m = 2 * fix( ( 2 * lmax + 1 ) / 2 ) + 1;  
du = 2 * pi / m;
u = du * ( 0 : m - 1 );
%  quadrature points and weights for polar direction
n = 2 * ( lmax + 2 );
[ x, w ] = lgwt( n, -1, 1 );
t = acos( x ); 

%  grid for quadrature points and weights
[ u, t ] = ndgrid( u, t );
w = du * ones( m, 1 ) * w .';
%  set output
if ~exist( 'key', 'var' ) || strcmp( key, 'matrix' )
  quad = struct( 'u', u, 't', t, 'w', w );
else
  quad = struct( 'u', u( : ), 't', t( : ), 'w', w( : ) );
end
