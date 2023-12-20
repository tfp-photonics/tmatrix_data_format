function sol = solve( obj, q, varargin )
%  Solve - Solve BEM equation using characteristic modes.
%
%  Usage for obj = charModes :
%    sol = solve( obj, q, PropertyPairs )
%  Input
%    q      :  structure with inhomogeneities and wavenumber
%  PropertyName
%    ind    :  index to selected modes
%  Output
%    sol    :  solution with tangential electromagnetic fields

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'ind', 1 : 2 * ndof( obj.bem.tau ) );
%  parse input
parse( p, varargin{ : } );

%  characteristic modes 
obj = eval( obj, q.k0 );
w = cat( 1, obj.vec.h, - 1i * obj.vec.e );
%  inhomogeneity for excitation
[ n, siz ] = deal( ndof( obj.bem.tau ), size( q.e ) );
q = reshape( cat( 1, q.e, 1i * q.h ), 2 * n, [] );

%  solution vector
a = bsxfun( @rdivide,  ...
  w( :, p.Results.ind ) .' * q, 1 + 1i * obj.val( p.Results.ind ) );
u = w( :, p.Results.ind ) * a;
%  set up solution object
sol = galerkin.solution( obj.bem.tau, obj.k0,  ...
  1i * reshape( u( n + 1 : end, : ), siz ), reshape( u( 1 : n, : ), siz ) );  
