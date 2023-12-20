function y = grid( obj, pts, varargin )
%  GRID - Tabulation grids for position pairs.
%
%  Usage for obj = stratified.layerstructure :
%    y = grid( obj, pts, PropertyPairs )
%  Input
%    pts    :  position pairs, see SLICE
%  PropertyName
%    nr     :  number of radial positions for tabulation
%    nz     :  number of z-positions for tabulation
%  Output
%    y      :  grid structure

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'nr', 20 );
addParameter( p, 'nz', 20 );
%  parse input
parse( p, varargin{ : } );

%  tabulation ranges for position pairs
y = range( obj, pts, varargin{ : } );
%  loop over ranges
for it = 1 : numel( y )
  y( it ).r  = linspace( y( it ).r(  1 ), y( it ).r(  2 ), p.Results.nr );
  y( it ).z1 = linspace( y( it ).z1( 1 ), y( it ).z1( 2 ), p.Results.nz );
  if ~isempty( y( it ).z2 )
    y( it ).z2 = linspace( y( it ).z2( 1 ), y( it ).z2( 2 ), p.Results.nz );
  end
end
