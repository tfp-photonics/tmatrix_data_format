function obj = rot( obj, angle )
%  ROT - Rotate polygon by given angle.
%
%  Usage for obj = polygon :
%    obj = rot( obj, angle )
%  Input
%    angle  :  rotation angle in degrees

%  transform angle to rad
angle = angle / 180 * pi;
%  rotation matrix
rot = [ cos( angle ), sin( angle ); - sin( angle ), cos( angle ) ];

%  rotate positions
for it = 1 : numel( obj )
  obj( it ).pos = obj( it ).pos * rot;
end
