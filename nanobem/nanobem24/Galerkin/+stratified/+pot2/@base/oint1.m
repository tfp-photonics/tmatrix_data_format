function data = oint1( obj, k0, varargin )
%  OINT1 - Default boundary integral evaluators.
%
%  Usage for obj = stratified.pot2.base :
%    data = oint1( obj, k0, PropertyPairs )
%  Input
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    ind    :  tensor indices [i1,i2,q2,a2,k]
%  Output
%    auxiliary information for evaluation of boundary integrals

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ind', 1 : 5 );
%  parse input
parse( p, varargin{ : } );

%  evaluation and integration points 
[ pt1, pt2 ] = deal( obj.pt1, obj.pt2 );
%  integration points, weights and shape functions
pos1 = pt1.pos;
[ pos2, w, f, fp ] = eval( pt2 );

%  save positions
data = struct( 'pos1', pos1, 'pos2', pos2 );
%  radial distance
pos1( :, 3 ) = 0;
pos2( :, :, 3 ) = 0;
r = pdist2( pos1, reshape( pos2, [], 3 ) );
r = reshape( r, size( pos1, 1 ), size( pos2, 1 ), [] );
%  avoid division by zero
r( r < 1e-10 ) = 1e-10;

%  dummy indices for tensor class
ind = num2cell( p.Results.ind );
[ i1, i2, q2, a2, k ] = deal( ind{ : } );
%  unit vectors
data.r = tensor( r, [ i1, i2, q2 ] );
data.er = ( tensor( pos1, [ i1, k ] ) - tensor( pos2, [ i2, q2, k ] ) ) ./ data.r;
data.ez = tensor( [ 0, 0, 1 ], k );
data.et = cross( data.ez, data.er, k ); 
%  shape function and derivative
data.f  = tensor( f,  [ i2, q2, a2, k ] ) * tensor( w, [ i2, q2 ] );
data.fp = tensor( fp, [ i2, q2, a2 ]    ) * tensor( w, [ i2, q2 ] );

%  material properties of layer materials
mat1 = obj.pt1.mat( obj.pt1.imat );
mat2 = obj.pt2.mat( obj.pt2.inout( 2 ) );
[ data.k1, data.mu1, data.eps1 ] = deal( mat1.k( k0 ), mat1.mu( k0 ), mat1.eps( k0 ) );
[ data.k2, data.mu2, data.eps2 ] = deal( mat2.k( k0 ), mat2.mu( k0 ), mat2.eps( k0 ) );
