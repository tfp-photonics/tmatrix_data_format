function [ t, w, ra, rb ] = adapt2( obj, t0, t1, t2, ta, da, tb, db )
%  ADAPT2 - Adapt integration points to triangles (aloof geometry).
%
%  Usage for obj = triquad : 
%    [ t, w, ra, rb ] = adapt2( obj, t0, t1, t2, ta, da, tb, db )
%  Input
%    t0       :  offset angle
%    t1,t2    :  range of angles
%    ta,da    :  parameters for first line
%    tb,db    :  parameters for second line
%  Output
%    t        :  integration angles
%    w        :  integration weight
%    ra       :  minimum radius 
%    rb       :  maximum radius

%  integration points and weights for unit range
[ xt, wt ] = lgwt( obj.npol( end ), 0, 1 );  

%  dummy indices for tensor class
[ i, qt ] = deal( 1, 2 );
%  convert input to tensor
[ xt, wt ] = deal( tensor( xt, qt ), tensor( wt, qt ) );
[ t1, t2 ] = deal( tensor( t1, i ), tensor( t2, i ) );
[ ta, da ] = deal( tensor( ta, i ), tensor( da, i ) );
[ tb, db ] = deal( tensor( tb, i ), tensor( db, i ) );

%  angles and radii
t = t1 + ( t2 - t1 ) * xt;
ra = da ./ cos( t - ta );
rb = db ./ cos( t - tb );
%  integration weights
w = ( t2 - t1 ) * wt;
t = t + tensor( t0, i );

%  convert integration points and weights to numeric
t = double( t, [ i, qt ] );
w = double( w, [ i, qt ] );
ra = double( ra, [ i, qt ] );
rb = double( rb, [ i, qt ] );
%  sort radii
i1 = ra > rb;
[ ra( i1 ), rb( i1 ) ] = deal( rb( i1 ), ra( i1 ) );
