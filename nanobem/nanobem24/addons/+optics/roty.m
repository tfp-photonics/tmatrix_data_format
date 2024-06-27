function rot = roty( t, key )
%  ROTY - Rotation around y-axis.
%
%  Usage :
%    rot = optics.roty( t, key )
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
rot = [ cost, 0, -sint; 0, 1, 0; sint, 0, cost ];
