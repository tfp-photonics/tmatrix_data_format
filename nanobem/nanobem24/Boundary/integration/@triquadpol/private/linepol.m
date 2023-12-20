function [ t0, d0 ] = linepol( t1, r1, t2, r2 )
%  LINEPOL - Line between vertices in polar coordinates.
%
%  Usage : 
%    [ t, r ] = linepol( t1, r1, t2, r2 )
%  Input
%    t1,r1  :  positions of first point in polar coordinates
%    t2,r2  :  positions of second point in polar coordinates
%  Output
%     t0,r0 :  r(t)=d0/cos(t-t0)

t2( t1 == t2 ) = t2( t1 == t2 ) + 1e-10;
sig = sign( sin( t2 - t1 ) );
R = sqrt( r1 .^ 2 + r2 .^ 2 - 2 .* r1 .* r2 .* cos( t2 - t1 ) );
%  line parameters
d0 = sig ./ R .* r1 .* r2 .* sin( t2 - t1 );
t0 = atan2( - sig .* ( r2 .* cos( t2 ) - r1 .* cos( t1 ) ),  ...
              sig .* ( r2 .* sin( t2 ) - r1 .* sin( t1 ) ) );
