function obj = secondary( obj, varargin )
%  SECONDARY - Secondary fields for layer structure.
%
%  Usage for obj = optics.decompose :
%    obj = secondary( obj, layer, PropertyPairs )
%  Input
%    layer    :  layer structure
%  PropertyName
%    dir      :  'up' for upgoing or 'down' for downgoing waves
%  Output
%    obj      :  plane wave decomposition of secondary fields

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'layer', [] );
addParameter( p, 'dir', 'up' );
%  parse input
parse( p, varargin{ : } );

layer = p.Results.layer;
%  upgoing and downgoing waves
i1 = obj.dir( :, 3 ) > 0;
i2 = obj.dir( :, 3 ) < 0;
%  set output
[ dir, efield ] = deal( [] );

%  upgoing waves
if nnz( i1 )
  %  reflected and transmitted waves 
  [ r, t ] = refl1( layer, obj.k0, obj.dir( i1, : ), obj.efield( i1, : ) );
  switch p.Results.dir
    case 'up'
      [ efield, dir ] = deal( [ efield; t.efield ], [ dir; t.dir ] );
    case 'down'
      [ efield, dir ] = deal( [ efield; r.efield ], [ dir; r.dir ] );
  end
end

%  downgoing waves
if nnz( i2 )
  %  reflected and transmitted waves 
  [ r, t ] = refl2( layer, obj.k0, obj.dir( i2, : ), obj.efield( i2, : ) );
  switch p.Results.dir
    case 'up'
      [ efield, dir ] = deal( [ efield; r.efield ], [ dir; r.dir ] );
    case 'down'
      [ efield, dir ] = deal( [ efield; t.efield ], [ dir; t.dir ] );
  end
end
%  set output
[ obj.efield, obj.dir ] = deal( efield, dir );
 