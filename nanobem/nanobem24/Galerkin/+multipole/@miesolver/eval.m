function t = eval( obj, k0 )
%  EVAL - Multipole coefficients for Mie theoty.
%
%  Usage for obj = multipole.tsolver :
%    t = eval( k0 )
%  Input
%    k0   :  wavenumber of light in vacuum
%  Output
%    t    :  T-matrix or multipole coefficients

%  Mie coefficients
[ a, b ] = miecoefficients( obj, k0, obj.tab.l );
%  initialize T-matrix
t = multipole.tmatrix( obj, k0 );
%  add matrices
n = numel( obj.tab.l );
t.aa = - diag( a );
t.ab =  zeros( n );
t.ba =  zeros( n );
t.bb = - diag( b );
