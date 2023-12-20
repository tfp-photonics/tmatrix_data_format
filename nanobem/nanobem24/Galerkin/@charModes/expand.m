function a = expand( obj, sol, varargin )
%  EXPAND - Expansion coefficients for characteristic modes.
%
%  Usage for obj = charModes :
%    a = expand( obj, sol, PropertyPairs )
%  Input
%    sol    :  BEM solution vector
%  PropertyName
%    ind    :  index to selected modes
%  Output
%    a      :  expansion coefficients

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'ind', 1 : 2 * ndof( sol.tau ) );
%  parse input
parse( p, varargin{ : } );

%  compute characteristic modes (if needed)
obj = eval( obj, sol.k0 );
%  solution vector and characteristic modes
u = cat( 1,     sol.h, - 1i *     sol.e );
w = cat( 1, obj.vec.h, - 1i * obj.vec.e );
%  reshape solution vector
siz = size( sol.e );
u = reshape( u, 2 * siz( 1 ), [] );

%  expansion coefficients
a = w( :, p.Results.ind ) .' * ( real( obj.A ) * u );
%  reshape output
a = reshape( a, [ size( a, 1 ), siz( 2 : end ) ] );
