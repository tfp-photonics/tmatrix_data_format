function [ e, h ] = fields( obj, pts, varargin )
%  FIELDS - Fields for planewave decomposition at requested points.
%
%  Usage for obj = optics.decompose :
%    [ e, h ] = fields( obj, pts, PropertyPairs )
%  Input
%    pts      :  points where electromagnetic fields are requested
%  PropertyName
%    imat     :  material index for unbounded medium
%    layer    :  layer structure
%    primary  :  reflected fields only
%  Output
%    e,h      :  electromagnetic fields

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'imat', 1 );
addParameter( p, 'layer', [] );
%  parse input
parse( p, varargin{ : } );

switch isempty( p.Results.layer )
  case 1
    [ e, h ] = fields1( obj, pts, p.Results.imat );
  otherwise
    [ e, h ] = fields2( obj, pts, p.Results.layer, varargin{ : } );
end
