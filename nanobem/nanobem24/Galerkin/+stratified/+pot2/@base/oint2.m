function data = oint2( ~, pt, varargin )
%  OINT2 - Refined boundary integral evaluators.
%
%  Usage for obj = stratified.pot2.base :
%    data = oint2( obj, pt, PropertyPairs )
%  Input
%    pt     :  polar integrator
%  PropertyName
%    ind    :  tensor indices [i,q,a,k]
%  Output
%    auxiliary information for evaluation of boundary integrals

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ind', 1 : 4 );
%  parse input
parse( p, varargin{ : } );

%  integration points, weights and shape functions
pos1 = pt.pos( pt.i1, : );
[ pos2, w, f, fp ] = eval( pt );
%  save positions
data = struct( 'pos1', pos1, 'pos2', pos2 );
pos1( :, 3 ) = 0;
pos2( :, :, 3 ) = 0;
%  radial distance
r = sqrt( bsxfun( @minus, pos1( :, 1 ), pos2( :, :, 1 ) ) .^ 2 +  ...
          bsxfun( @minus, pos1( :, 2 ), pos2( :, :, 2 ) ) .^ 2 );
%  avoid division by zero
r( r < 1e-10 ) = 1e-10;

%  dummy indices for tensor class
ind = num2cell( p.Results.ind );
[ i, q, a, k ] = deal( ind{ : } );
%  unit vectors
data.r = tensor( r, [ i, q ] );
data.er = ( tensor( pos1, [ i, k ] ) - tensor( pos2, [ i, q, k ] ) ) ./ data.r;
data.ez = tensor( [ 0, 0, 1 ], k );
data.et = cross( data.ez, data.er, k ); 
%  shape function and derivative
data.f  = tensor( f,  [ i, q, a, k ] ) * tensor( w, [ i, q ] );
data.fp = tensor( fp, [ i, q, a ]    ) * tensor( w, [ i, q ] );
