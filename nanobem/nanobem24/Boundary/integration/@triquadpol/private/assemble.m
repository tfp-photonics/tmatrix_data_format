function obj = assemble( obj, verts, pos, i1, t, w, ra, rb )
%  ASSEMBLE - Assemble quadrature rules for polar triangle integration.
%
%  Usage for obj = triquad :
%    quad = rquad( obj, verts, orgin, i1, t, w, ra, rb )
%  Input
%    verts    :  triangular vertices
%    pos      :  origins for polar integration
%    i1       :  indices for origins
%    t        :  discretized angles
%    w        :  integration weights
%    ra,rb    :  limits for radial integration

%  dummy indices for internal tensor class
[ i, qr, qt ] = deal( 1, 2, 3 );
%  convert to internal tensor class
t  = tensor( t,  [ i, qt ] );
w  = tensor( w,  [ i, qt ] );
ra = tensor( ra, [ i, qt ] );
rb = tensor( rb, [ i, qt ] );
%  integration points along radial direction
[ xr, wr ] = lgwt( obj.npol( 1 ), 0, 1 );
[ xr, wr ] = deal( tensor( xr, qr ), tensor( wr, qr ) );

%  expand integration points and weights along radial direction
r = ra + xr * ( rb - ra );
w =  w * wr * ( rb - ra ) * r;
t =  t * tensor( ones( 1, obj.npol( 1 ) ), qr );
%  origin positions
xi = tensor( pos( i1, 1 ), i );
yi = tensor( pos( i1, 2 ), i );
%  convert integration points to Cartesian coordinates
xr = r * cos( t ) + xi;
yr = r * sin( t ) + yi; 

%  triangle vertices
vx = verts( i1, :, 1 );
vy = verts( i1, :, 2 );
%  triangle corners
[ x1, y1 ] = deal( tensor( vx( :, 1 ), i ), tensor( vy( :, 1 ), i ) );
[ x2, y2 ] = deal( tensor( vx( :, 2 ), i ), tensor( vy( :, 2 ), i ) );
[ x3, y3 ] = deal( tensor( vx( :, 3 ), i ), tensor( vy( :, 3 ), i ) );
%  convert from Cartesian to triangle coordinates
a = ( y2 - y3 ) * ( x1 - x3 ) + ( x3 - x2 ) * ( y1 - y3 );
xt = ( ( y2 - y3 ) * ( xr - x3 ) + ( x3 - x2 ) * ( yr - y3 ) ) ./ a;
yt = ( ( y3 - y1 ) * ( xr - x3 ) + ( x1 - x3 ) * ( yr - y3 ) ) ./ a;

%  quadrature rules
x = reshape( double( xt, [ i, qt, qr ] ), numel( i1 ), [] );
y = reshape( double( yt, [ i, qt, qr ] ), numel( i1 ), [] );
w = reshape( double( w,  [ i, qt, qr ] ), numel( i1 ), [] );
%  set output
obj.quad = struct( 'x', x, 'y', y, 'w', w );
obj.i1 = i1;
