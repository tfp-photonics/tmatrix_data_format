function obj = mtimes( obj, rot )
%  MTIMES - Multiply with rotation matrix.
%
%  Usage for obj = optics.decompose :
%    obj = obj * rot
%    obj = rot * obj
%  Input
%    rot    :  rotation matrix
%  Output
%    obj    :  planewave decomposition with rotated fields and directions

if isa( rot, 'optics.decompose' )
  [ obj, rot ] = deal( rot, obj .' );
end
%  rotate fields and propagation directions
obj.efield = obj.efield * rot;
obj.dir = obj.dir * rot;
