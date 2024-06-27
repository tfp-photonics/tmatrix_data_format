function obj = times( obj, x )
%  TIMES - Multiply electric field with vector.
%
%  Usage for obj = optics.decompose :
%    obj = obj .* x
%  Input
%    x      :  vector
%  Output
%    obj    :  planewave decomposition multipled with vector

if isa( x, 'optics.decompose' )
  [ obj, x ] = deal( x, obj );
end
%  multiply field with vector
obj.efield = obj.efield .* x;
