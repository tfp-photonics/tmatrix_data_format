function [ e, x ] = efield( obj, field, varargin )
%  EFIELD - Electric field on image side.
%
%  Usage for obj = optics.lensimage2 :
%    [ e, x ] = efield( obj, field, PropertyPairs )
%  Input
%    field  :  electric far-fields or planewave decomposition
%  PropertyName
%    n      :  size of output array
%    focus  :  focus position of imaging lens
%  Output
%     e     :  electric image fields in focal plane
%     x     :  image coordinates

switch class( field )
  case 'double'
    [ e, x ] = efield1( obj, field, varargin{ : } );
  case 'optics.decompose'
    [ e, x ] = efield2( obj, field, varargin{ : } );
end
