function y = unitproject( x, u1, u2 )
%  UNITPROJECT - Project vector on unit vector
%
%  Usage :
%    y = unitproject( x, u1, u2 )
%  Input
%    x    :  vector
%    u1   :  unit vector
%    u2   :  second unit vector (optional)
%  Output
%    y    :  x.u1 or (x.u1)*u2

%  inner product
y = squeeze( sum( bsxfun( @times, x, u1 ), 2 ) );
%  expand with second unit vector?
if exist( 'u2', 'var' )
  y = double( tensor( y, [ 1, 2 ] ) * tensor( u2, [ 1, 3 ] ), [ 1, 3, 2 ] );
end
