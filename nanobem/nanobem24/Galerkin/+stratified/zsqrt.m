function y = zsqrt( x )
%  ZSQRT - Square root for complex argument.
%  Chose sign such that imaginary part is always positive.

y = sqrt( x );  
y = y .* sign( imag( y + 1i * eps ) );
