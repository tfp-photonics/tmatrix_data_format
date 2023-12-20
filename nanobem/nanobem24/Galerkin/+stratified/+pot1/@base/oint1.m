function [ oint, data ] = oint1( obj, k0 )
%  OINT1 - Defaults boundary integral evaluators.
%
%  Usage for obj = stratified.pot1.base :
%    [ oint, data ] = oint1( obj, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum
%  Output
%    oint   :  boundary integral evaluators
%    data   :  positions and material properties

%  integration points 
[ pt1, pt2 ] = deal( obj.pt1, obj.pt2 );
%  integration points, weights and shape functions
[ pos1, w1, f1, fp1 ] = eval( pt1 );
[ pos2, w2, f2, fp2 ] = eval( pt2 );
%  layer coordinates and distance
r = pdist2(  ...
  reshape( pos1( :, :, 1 : 2 ), [], 2 ), reshape( pos2( :, :, 1 : 2 ), [], 2 ) );
r( r < 1e-10 ) = 1e-10;
%  auxiliary data
data1 = struct( 'r', r, 'x1', pos1( :, :, 1 ), 'y1', pos1( :, :, 2 ),  ...
                        'x2', pos2( :, :, 1 ), 'y2', pos2( :, :, 2 ) );  

%  dummy indices for tensor class
[ i, q, a, k ] = deal( 1, 2, 3, 4 );
%  multiply shape functions with integration weights
f1 = double( tensor( f1, [ i, q, a, k ] ) * tensor( w1, [ i, q ] ), [ i, q, a, k ] );
f2 = double( tensor( f2, [ i, q, a, k ] ) * tensor( w2, [ i, q ] ), [ i, q, a, k ] );
fp1 = double( tensor( fp1, [ i, q, a ] ) * tensor( w1, [ i, q ] ), [ i, q, a ] );
fp2 = double( tensor( fp2, [ i, q, a ] ) * tensor( w2, [ i, q ] ), [ i, q, a ] );
%  Cartesian components
fx1 = f1( :, :, :, 1 );  fx2 = f2( :, :, :, 1 );
fy1 = f1( :, :, :, 2 );  fy2 = f2( :, :, :, 2 );
fz1 = f1( :, :, :, 3 );  fz2 = f2( :, :, :, 3 );

%  integration function and outer product
ifun  = @( f, g ) sum( reshape( f, size( g ) ) .* g, [ 2, 4 ] );
outer = @( f, g ) reshape( f, [], 1 ) * reshape( g, 1, [] );
%  boundary integral evaluators for single layer potential
oint.pp = @( a1, a2, g ) ifun( outer( fp1( :, :, a1 ), fp2( :, :, a2 ) ), g );
oint.zp = @( a1, a2, g ) ifun( outer( fz1( :, :, a1 ), fp2( :, :, a2 ) ), g );
oint.pz = @( a1, a2, g ) ifun( outer( fp1( :, :, a1 ), fz2( :, :, a2 ) ), g );
oint.zz = @( a1, a2, g ) ifun( outer( fz1( :, :, a1 ), fz2( :, :, a2 ) ), g );
oint.ss = @( a1, a2, g ) ifun( outer( fx1( :, :, a1 ), fx2( :, :, a2 ) ) +  ...
                               outer( fy1( :, :, a1 ), fy2( :, :, a2 ) ), g );
%  boundary integral evaluators for double layer potential
oint.rp = @( a1, a2, g ) ifun( triple1( a1, a2, fx1, fy1, fp2, data1 ), g ); 
oint.pr = @( a1, a2, g ) ifun( triple2( a1, a2, fp1, fx2, fy2, data1 ), g ); 
oint.rz = @( a1, a2, g ) ifun( triple1( a1, a2, fx1, fy1, fz2, data1 ), g ); 
oint.zr = @( a1, a2, g ) ifun( triple2( a1, a2, fz1, fx2, fy2, data1 ), g ); 

%  auxiliary data for evaluation of boundary integrals
data = struct( 'pos1', pos1, 'pos2', pos2 );
%  material properties of layer material
mat1 = obj.pt1.mat( obj.pt1.inout( 2 ) );
mat2 = obj.pt2.mat( obj.pt2.inout( 2 ) );
[ data.k1, data.mu1, data.eps1 ] = deal( mat1.k( k0 ), mat1.mu( k0 ), mat1.eps( k0 ) );
[ data.k2, data.mu2, data.eps2 ] = deal( mat2.k( k0 ), mat2.mu( k0 ), mat2.eps( k0 ) );


function val = triple1( a1, a2, fx, fy, g, data1 )
%  TRIPLE1 - [ fx * ( y1 - y2 ) - fy * ( x1 - x2 ) ] * g.
outer = @( f, g ) reshape( f, [], 1 ) * reshape( g, 1, [] );
val = ( outer( fx( :, :, a1 ) .* data1.y1 -  ...
               fy( :, :, a1 ) .* data1.x1, g( :, :, a2 ) ) -  ...
        outer( fx( :, :, a1 ), g( :, :, a2 ) .* data1.y2 ) +  ...
        outer( fy( :, :, a1 ), g( :, :, a2 ) .* data1.x2 ) ) ./ data1.r;


function val = triple2( a1, a2, g, fx, fy, data1 )
%  TRIPLE2 - g * [ fx * ( y1 - y2 ) - fy * ( x1 - x2 ) ].
outer = @( g, f ) reshape( g, [], 1 ) * reshape( f, 1, [] );
val = ( outer( g( :, :, a1 ) .* data1.y1, fx( :, :, a2 ) ) -  ...
        outer( g( :, :, a1 ) .* data1.x1, fy( :, :, a2 ) ) -  ...
        outer( g( :, :, a1 ), fx( :, :, a2 ) .* data1.y2   -  ...
                              fy( :, :, a2 ) .* data1.x2 ) ) ./ data1.r;
                            