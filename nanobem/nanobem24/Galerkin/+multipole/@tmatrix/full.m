function t = full( obj )
%  FULL - Expand T-matrices to full size.
%
%  Usage for obj = multipole.tmatrix :
%    t = full( obj )
%  Output
%    t    :  full T-matrices

switch numel( obj )
  case 1
      t = [ obj.aa, obj.ab; obj.ba, obj.bb ];
  otherwise
    t = arrayfun( @( x ) full( x ), obj, 'uniform', 0 );
    t = cat( 3, t{ : } );
end
