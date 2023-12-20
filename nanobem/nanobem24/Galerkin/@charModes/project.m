function [ sol, a ] = project( obj, sol, varargin )
%  PROJECT - Project BEM solution on characteristic modes.
%
%  Usage for obj = charModes :
%    [ sol, a ] = project( obj, sol, PropertyPairs )
%  Input
%    sol    :  BEM solution vector
%  PropertyName
%    ind    :  index to selected modes
%  Output
%    sol    :  BEM solution projected on modes
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

%  expansion coefficient and mode projection
a = w( :, p.Results.ind ) .' * ( real( obj.A ) * u );
u = w( :, p.Results.ind ) * a;
%  set output vector
sol.e = 1i * reshape( u( 1 + siz( 1 ) : end, : ), siz );
sol.h =      reshape( u( 1 : siz( 1 ),       : ), siz );
