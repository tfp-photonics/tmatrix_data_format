function y = laguerre( n, m, x )
%  LAGUERRE - Associated Laguerre polynomial.

switch n
  case 0
    y = 1;
  otherwise
    %  compute polynomials
    P = zeros( n + 1, 1 );
    for u = 0 : n
      P( n + 1 - u ) = ( -1 ) ^ u * factorial( n + m ) / ...
        ( factorial( n - u ) * factorial( m + u ) * factorial( u ) );
    end
    %  evaluate polynomials
    y = polyval( P, x );
end
