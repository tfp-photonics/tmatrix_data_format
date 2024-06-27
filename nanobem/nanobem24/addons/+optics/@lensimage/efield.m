function e = efield( obj, field, varargin )
%  EFIELD - Electric field on image side, Eq. (3.10).
%
%  Usage for obj = optics.lensimage :
%    e = efield( obj, field, x, y, PropertyPairs )
%  Input
%    field  :  electric far-fields or planewave decomposition
%    x,y    :  image coordinates
%  PropertyName
%    focus  :  focus position of imaging lens
%  Output
%     e     :  electric image fields in focal plane

switch class( field )
  case 'double'
    e = efield1( obj, field, varargin{ : } );
  case 'optics.decompose'
    e = efield2( obj, field, varargin{ : } );
end
