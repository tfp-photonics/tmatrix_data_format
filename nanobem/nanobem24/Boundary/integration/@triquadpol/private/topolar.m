function [ t, t0, r ] = topolar( verts, orgin )
%  TOPOLAR - Triangle vertices in polar coordinates.
%
%  Usage : 
%    [ t, t0, r ] = topolar( verts, orgin )
%  Input
%    verts  :  triangle vertices
%    orgin  :  origins for polar integration
%  Output
%    t      :  angles of triangle vertices
%    t0     :  offset angle
%    r      :  radii of triangle vertices

%  angle between impact positions and boundary centroids
pos = reshape( squeeze( sum( verts, 2 ) ) / 3, [], 2 );
t0 = cart2pol( pos( :, 1 ) - orgin( :, 1 ), pos( :, 2 ) - orgin( :, 2 ) );
%  expand vertices in rotated system
xx = bsxfun( @minus, verts( :, :, 1 ), orgin( :, 1 ) );
yy = bsxfun( @minus, verts( :, :, 2 ), orgin( :, 2 ) );
x =   bsxfun( @times, xx, cos( t0 ) ) + bsxfun( @times, yy, sin( t0 ) );
y = - bsxfun( @times, xx, sin( t0 ) ) + bsxfun( @times, yy, cos( t0 ) );
%  convert to polar coordinates
[ t, r ] = cart2pol( x, y );

%  sort according to angle
[ t, i1 ] = sort( t, 2 );
for it = 1 : size( verts, 1 )
  r( it, : ) = r( it, i1( it, : ) );
end
