function rot = rotx( t, key )
%  ROTX - Rotation around x-axis.
%
%  Usage :
%    rot = optics.rotx( t, key )
%  Input
%    t    :  rotation angle
%    key  :  'deg' or 'rad'
%  Output
%    rot  :  rotation matrix

%  rotation angle
if ~exist( 'key', 'var' ) || strcmp( key, 'deg' )
  t = t / 180 * pi;
end
%  rotation matrix
[ sint, cost ] = deal( sin( t ), cos( t ) );
rot = [ 1, 0, 0; 0, cost, -sint; 0, sint, cost ];
