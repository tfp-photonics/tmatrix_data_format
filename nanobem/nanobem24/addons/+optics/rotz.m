function rot = rotz( t, key )
%  ROTZ - Rotation around z-axis.
%
%  Usage :
%    rot = optics.rotz( t, key )
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
rot = [ cost, -sint, 0; sint, cost, 0; 0, 0, 1 ];
