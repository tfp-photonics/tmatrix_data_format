function obj = toggle( obj, inout )
%  TOGGLE - Toggle between inside and outside fields.
%
%  Usage for obj = feibelman.solution :
%    obj = toggle( obj, inout )
%  Input
%    inout  :  store in E,H fields at particle inside or outside

%  toggle fields ?
if obj.inout ~= inout
  [ obj.e, obj.et ] = deal( obj.et, obj.e );
  [ obj.h, obj.ht ] = deal( obj.ht, obj.h );
end
%  set INOUT flag of solution
obj.inout = inout;
