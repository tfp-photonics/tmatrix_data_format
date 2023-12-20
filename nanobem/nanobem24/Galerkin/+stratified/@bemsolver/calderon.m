function cal = calderon( obj, k0 )
%  CALDERON - Calderon matrix.
%
%  Usage for stratified.bemsolver :
%    cal = calderon( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum
%  Output
%    cal  :  Calderon matrix

%  fill tabulated Green functions
obj.green = fill( obj.green, k0 );
%  compute direct and reflected Calderon matrices
siz = ndof( obj.tau ) * [ 1, 1 ];
cal = calderon( obj.pot1, k0 ) +  ...
      calderon( obj.pot2, obj.green, k0, 'siz', siz );
    