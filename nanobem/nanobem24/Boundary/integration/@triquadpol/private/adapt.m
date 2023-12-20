function obj = adapt( obj )
%  ADAPT - Adapt polar integration points for triangles.

%  rotate boundary elements and origins to triangle planes
[ verts, pos ] = to2d( obj );
%  centroid position and bounding box radius
pos2 = squeeze( sum( verts, 2 ) ) / 3;
rad = boxradius( obj.tau );
% origins located inside of bounding box circle ?
u = pos2 - pos;
in = sqrt( dot( u, u, 2 ) ) <= rad( : );
%  loop over tagged pairs
for i1 = reshape( find( in ), 1, [] )
  in( i1 ) = inpolygon(  ...
    pos( i1, 1 ), pos( i1, 2 ), verts( i1, :, 1 ), verts( i1, :, 2 ) );
end

%  origins inside of triangles
if nnz( in )
  %  convert triangle vertices to polar coordinates
  [ t, t0, r ] = topolar( verts( in, :, : ), pos( in, :, : ) );
  %  line in polar coordinates
  [ t12, d12 ] = linepol( t( :, 1 ), r( :, 1 ), t( :, 2 ), r( :, 2 ) );
  [ t23, d23 ] = linepol( t( :, 2 ), r( :, 2 ), t( :, 3 ), r( :, 3 ) );
  [ t31, d31 ] = linepol( t( :, 3 ), r( :, 3 ), t( :, 1 ), r( :, 1 ) );  
  %  adapt integration points to triangle
  [ t1, w1, ra1, rb1 ] = adapt1( obj, t0, t( :, 1 ), t( :, 2 ), t12, d12 );
  [ t2, w2, ra2, rb2 ] = adapt1( obj, t0, t( :, 2 ), t( :, 3 ), t23, d23 );  
  [ t3, w3, ra3, rb3 ] = adapt1( obj, t0, t( :, 3 ), t( :, 1 ), t31, d31 ); 
  %  assemble integration points  
  yout{ 1 } = assemble( obj, verts, pos, find( in ),  ...
    [ t1, t2, t3 ], [ w1, w2, w3 ], [ ra1, ra2, ra3 ], [ rb1, rb2, rb3 ] );
end

%  origins outside of triangles
if nnz( ~in )
  %  convert triangle vertices to polar coordinates
  [ t, t0, r ] = topolar( verts( ~in, :, : ), pos( ~in, :, : ) );
  %  line in polar coordinates
  [ t12, d12 ] = linepol( t( :, 1 ), r( :, 1 ), t( :, 2 ), r( :, 2 ) );
  [ t13, d13 ] = linepol( t( :, 1 ), r( :, 1 ), t( :, 3 ), r( :, 3 ) );
  [ t23, d23 ] = linepol( t( :, 2 ), r( :, 2 ), t( :, 3 ), r( :, 3 ) );
  %  adapt integration points to triangle
  [ t1, w1, ra1, rb1 ] = adapt2( obj, t0, t( :, 1 ), t( :, 2 ), t12, d12, t13, d13 );
  [ t2, w2, ra2, rb2 ] = adapt2( obj, t0, t( :, 2 ), t( :, 3 ), t23, d23, t13, d13 );
  %  assemble integration points  
  yout{ 2 } = assemble( obj, verts, pos, find( ~in ),  ...
    [ t1, t2 ], [ w1, w2 ], [ ra1, ra2 ], [ rb1, rb2 ] );
end

%  set output
obj = horzcat( yout{ : } );
