function obj = fill( obj, k0 )
%  FILL - Fill matrices and compute Calderon matrix.
%
%  Usage for stratified.bemsolver :
%    obj = fill( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum

%  compute Calderon matrix only if needed
if isempty( obj.k0 ) || obj.k0 ~= k0
  %  fill tabulated Green functions
  obj.green = fill( obj.green, k0 );
  %  compute direct and reflected Calderon matrices
  siz = ndof( obj.tau ) * [ 1, 1 ];
  cal = calderon( obj.pot1, k0 ) +  ...
        calderon( obj.pot2, obj.green, k0, 'siz', siz );
  %  LU decomposition of Calderon matrix
  obj.cal = decomposition( cal, 'LU' );
  obj.k0 = k0;
end
