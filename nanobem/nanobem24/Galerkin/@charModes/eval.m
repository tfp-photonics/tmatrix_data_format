function obj = eval( obj, k0 )
%  EVAL - Evaluate characteristic modes.
%
%  Usage for obj = charModes :
%    obj = eval( obj, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum

if isempty( obj.k0 ) || obj.k0 ~= k0
  %  Calderon matrix
  obj.k0 = k0;
  A = calderon( obj.bem, k0 );
  %  subdivide matrix
  [ tau, n ] = deal( obj.bem.tau, ndof( obj.bem.tau ) );
  A = mat2cell( A, [ n, n ], [ n, n ] );
  %  Calderon matrix for characteristic modes
  A = [ A{ 1, 2 }, 1i * A{ 1, 1 }; 1i * A{ 2, 2 }, - A{ 2, 1 } ];
  
  %  symmetrized Calderon matrix and symmetry number
  obj.A = 0.5 * ( A + A .' );
  obj.rsym = norm( A - A .', 'fro' ) / norm( obj.A, 'fro' );
  %  eigenvalues and characteristic modes
  [ w, val ] = eig( imag( obj.A ), real( obj.A ), 'vector' );
  %  sort according to eigenvalue
  [ ~, ind ] = sort( val );
  [ obj.val, w ] = deal( val( ind ), w( :, ind ) );
  %  normalization
  w = bsxfun( @rdivide, w, sqrt( diag( w .' * real( obj.A ) * w ) ) .' );
  
  %  set up eigenvectors
  obj.vec = galerkin.solution(  ...
    tau, k0, 1i * w( n + 1 : end, : ), w( 1 : n, : ) );  
end
