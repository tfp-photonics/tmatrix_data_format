function obj = transform( obj, varargin )
%  TRANSFORM - Scale, shift, or flip particle boundary.
%
%  Usage for obj = particle :
%    obj = plot( obj, PropertyPair )
%  PropertyName
%    'scale'    :  scale  particle boundary
%    'shift'    :  shift  particle boundary
%    'rot'      :  rotate particle boundary
%    'flip'     :  flip   particle boundary
%  Output
%    obj        :  transformed particle boundary

for i = 1 : 2 : numel( varargin )
  val = varargin{ i + 1 };
  switch varargin{ i }
    case 'scale'
      obj.verts = obj.verts .* val;
    case 'shift'
      obj.verts = obj.verts + val;
    case 'rot'
      obj.verts = obj.verts * val .';
    case 'flip'
      obj.verts( :, val ) = - obj.verts( :, val );
      obj.faces = fliplr( obj.faces );
  end
end