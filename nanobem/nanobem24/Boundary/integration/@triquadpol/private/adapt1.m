function [ t, w, ra, rb ] = adapt1( obj, t0, t1, t2, tb, db )
%  ADAPT1 - Adapt integration points to triangles (penetrating beams).
%
%  Usage for obj = triquadpol : 
%    [ t, w, ra, rb ] = adapt1( obj, t0, t1, t2, tb, db )
%  Input
%    t0       :  offset angle
%    t1,t2    :  range of angles
%    tb,db    :  line parameters
%  Output
%    t        :  integration angles
%    w        :  integration weight
%    ra       :  minimum radius 
%    rb       :  maximum radius

%  integration points and weights for unit range
[ xt, wt ] = lgwt( obj.npol( end ), 0, 1 );  
%  assure proper range of angles
t2( t2 < t1 ) = t2( t2 < t1 ) + 2 * pi;

%  dummy indices for tensor class
[ i, qt ] = deal( 1, 2 );
%  convert input to tensor
[ xt, wt ] = deal( tensor( xt, qt ), tensor( wt, qt ) );
[ t1, t2 ] = deal( tensor( t1, i ), tensor( t2, i ) );
[ tb, db ] = deal( tensor( tb, i ), tensor( db, i ) );

%  angles and radii
t = t1 + ( t2 - t1 ) * xt;
ra = 0 * t;
rb = db ./ cos( t - tb );
%  integration weights
w = ( t2 - t1 ) * wt;
t = t + tensor( t0, i );

%  convert integration points and weights to numeric
t = double( t, [ i, qt ] );
w = double( w, [ i, qt ] );
ra = double( ra, [ i, qt ] );
rb = double( rb, [ i, qt ] );
