function [ up, ut, ur, p, t ] = cart2unit( vec )
%  CART2UNIT - Unit vectors.
%
%  Usage :
%    [ up, ut, ur, p, t ] = cart2unit( vec )
%  Input
%    vec    :  vectors
%  Output
%    up     :  unit vector in azimuthal direction  (phi)
%    ut     :  unit vector in polar direction      (theta)
%    ur     :  unit vector in radial direction     (rho)
%    p      :  azimuthal angle
%    t      :  polar angle

siz = size( vec );
vec = reshape( vec, [], 3 );
%  azimuthal and polar angles
[ p, t ] = cart2sph( vec( :, 1 ), vec( :, 2 ), vec( :, 3 ) );
t = pi / 2 - t;

[ sinp, cosp ] = deal( sin( p ), cos( p ) );
[ sint, cost ] = deal( sin( t ), cos( t ) );
% unit vectors
up = [ - sinp, cosp, 0 * p ];
ur = [   cosp, sinp, 0 * p ];
ut = [   cosp .* cost, sinp .* cost, - sint ];

%  reshape output
up = reshape( up, siz );
ur = reshape( ur, siz );
ut = reshape( ut, siz );
