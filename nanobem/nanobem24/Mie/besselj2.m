function  j = besselj2( n, x )
%  BESSELJ2 - Fast Bessel function using recursion.
%
%  Usage :
%    j = besselj2( n, x )
%  Input
%    n    :  array of equation orders
%    x    :  arguments for Bessel functions
%  Output
%    j    :  Bessel functions

%  reshape argument
x = reshape( x, 1, [] );
%  allocate table of Bessel functions
nmax = max( max( abs( n ) ), 1 );
ntab = 0 : nmax;
jtab = zeros( nmax + 1, numel( x ) );
%  downward recursion
jtab( nmax + 1, : ) = besselj( nmax,     x );
jtab( nmax,     : ) = besselj( nmax - 1, x );
%  loop over evaluation orders
for n1 = nmax - 1 : -1 : 1
  jtab( n1, : ) = 2 * n1 ./ x .* jtab( n1 + 1, : ) - jtab( n1 + 2, : );
end
%  deal with zeros
ind = abs( x ) < 1e-10;
jtab( :, ind ) = 0;
jtab( 1, ind ) = besselj( 0, x( ind ) );

%  allocate output
j = zeros( numel( n ), numel( x ) );
%  set output for positive orders
[ ind, i1 ] = ismember( + n, ntab );  
j( ind, : ) = jtab( i1( ind ), : );
%  set output for negative orders
[ ind, i2 ] = ismember( - n, ntab);  
if any( ind )
  j( ind, : ) = bsxfun( @times, jtab( i2( ind ), : ), ( - 1 ) .^ ntab( i2( ind ) ) .' );
end
