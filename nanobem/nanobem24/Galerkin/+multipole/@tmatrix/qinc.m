function q = qinc( obj, fun, varargin )
%  QINC - Mie coefficients for incoming excitation.
%
%  Usage for obj = multipole.tmatrix :
%    q = qinc( obj, fun, PropertyPairs )
%  Input
%    fun        :  incoming fields [e,h]=fun(pos,k0) 
%  PropertyName
%    diameter   :  diameter for multipole expansion
%  Output
%    q          :  structure with incoming multipole coefficients

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'diameter', 100 );
%  parse input
parse( p, varargin{ : } );

%  table of angular degrees and orders, wavenumber of light in vacuum
[ tab, k0 ] = deal( obj.solver.tab, obj.k0 );
%  quadrature points for azimuthal direction
lmax = max( tab.l );
n1 = 2 * fix( ( 2 * lmax + 1 ) / 2 ) + 1;
x1 = ( 0 : n1 - 1 ) / n1 * 2 * pi;
%  quadrature points and weights for polar direction
n2 = 2 * ( lmax + 2 );
[ x2, w2 ] = lgwt( n2, 0, pi );

%  sphere radius and outer product
r = 0.5 * p.Results.diameter;
outer = @( y1, y2 ) reshape( y1, [], 1 ) * reshape( y2, 1, [] );
%  integration positions
pos = cat( 3, r * outer( cos( x1 ), sin( x2 ) ),  ...
              r * outer( sin( x1 ), sin( x2 ) ),  ...
              r * outer(   x1 .^ 0, cos( x2 ) ) );
%  incoming electromagnetic fields
[ e, h ] = fun( reshape( pos, [], 3 ), k0 );
n3 = numel( e ) / numel( pos );
%  dot product with position vector
e = bsxfun( @times, reshape( e, [], n3 ), pos( : ) );
h = bsxfun( @times, reshape( h, [], n3 ), pos( : ) );
e = reshape( sum( reshape( e, [], 3, n3 ), 2 ), n1, [] );
h = reshape( sum( reshape( h, [], 3, n3 ), 2 ), n1, [] );

%  FFT for azimuthal angle
e = fftshift( fft( e, [], 1 ), 1 ) * 2 * pi / n1;
h = fftshift( fft( h, [], 1 ), 1 ) * 2 * pi / n1;
%  expand to angular degrees
[ ~, ind ] = ismember( tab.m, - fix( 0.5 * n1 ) : fix( 0.5 * n1 ) );
e = reshape( e( ind, : ), [], n3 );
h = reshape( h( ind, : ), [], n3 );
%  Legendre polynomials and prefactors
y = spharm( tab.l, tab.m, x2, 0 * x2 );
y = bsxfun( @times, y, reshape( w2 .* sin( x2 ), 1, [] ) );
y = bsxfun( @rdivide, y, sqrt( tab.l( : ) .* ( tab.l( : ) + 1 ) ) );
%  integration with Legendre polynomials
e = sum( reshape( bsxfun( @times, e, y( : ) ), [], n2, n3 ), 2 );
h = sum( reshape( bsxfun( @times, h, y( : ) ), [], n2, n3 ), 2 );

%  material parameters and Bessel function
mat = obj.solver.embedding;
[ k, Z ] = deal( mat.k( k0 ), mat.Z( k0 ) );
j = riccatibessel( tab.l, r * k, 'j' );
%  set output
q.a = - bsxfun( @rdivide, squeeze( e ), j( : ) ) * k / Z;
q.b =   bsxfun( @rdivide, squeeze( h ), j( : ) ) * k;
%  set wavenumber
q.k0 = k0;
