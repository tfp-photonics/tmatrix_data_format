function obj = fill( obj, k0 )
%  FILL - Fill Calderon matrix.
%
%  Usage for obj = galerkin.bemsolver :
%    obj = fill( obj, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum
%  Output
%    obj    :  BEM solver with pre-computed Calderon matrix

if isempty( obj.cal ) || obj.k0 ~= k0
  %  LU decomposition of Calderon matrix
  cal = calderon( obj.pot, k0 );
  obj.cal = decomposition( cal, 'lu' );
  %  save wavenumber
  obj.k0 = k0;
end
