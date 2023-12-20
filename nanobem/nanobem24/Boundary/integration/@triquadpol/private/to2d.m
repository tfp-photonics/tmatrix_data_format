function [ verts, pos ] = to2d( obj )
%  TO2D - Rotate to triangle plane.
%
%  Usage : 
%    [ verts, pos ] = to2d( obj )
%  Output
%    verts  :  rotated 2d coordinates
%    pos    :  rotated origin positions

%  extract input
[ tau, pos ] = deal( obj.tau, obj.pos );
%  triangle vertices
verts = permute( cat( 3, tau.verts ), [ 3, 1, 2 ] );
v1 = squeeze( verts( :, 1, : ) );
v2 = squeeze( verts( :, 2, : ) );
v3 = squeeze( verts( :, 3, : ) );
%  difference vectors
v31 = v3 - v1;
v21 = v2 - v1;

%  normalization function
norm = @( x ) bsxfun( @rdivide, x, sqrt( dot( x, x, 2 ) ) );
%  unit vectors
r1 = norm( v21 );
r3 = norm( cross( r1, v31, 2 ) );
r2 = cross( r3, r1, 2 );

%  vertices in triangle plane
v2 = [ dot( r1, v21, 2 ), dot( r2, v21, 2 ) ];
v3 = [ dot( r1, v31, 2 ), dot( r2, v31, 2 ) ];
verts = permute( cat( 3, 0 * v2, v2, v3 ), [ 1, 3, 2 ] );
%  origin in triangle plane
pos = cat( 2, dot( r1, pos - v1, 2 ), dot( r2, pos - v1, 2 ) );
