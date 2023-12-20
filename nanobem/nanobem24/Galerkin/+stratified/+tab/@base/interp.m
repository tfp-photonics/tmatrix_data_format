function y = interp( obj, pos1, pos2, k0, varargin )
%  INTERP - Interpolate Green function table.
%
%  Usage for obj = stratified.tab.base :
%    y = interp( obj, pos1, pos2, k0, PropertyPairs )
%  Input
%    pos1     :  observation points, must be all located in same layer
%    pos2     :  source points, must be all located in same layer
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    rows     :  positions of same size or not
%    smooth   :  w/o or with quasistatic contribution
%    stat     :  quasistatic contribution only
%  Output
%    y        :  Green function table

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'stat', 0 );
%  parse input
parse( p, varargin{ : } );

%  layer indices
i1 = unique( indlayer( obj( 1 ).layer, reshape( pos1, [], 3 ) ) );
i2 = unique( indlayer( obj( 1 ).layer, reshape( pos2, [], 3 ) ) );
assert( isscalar( i1 ) && isscalar( i2 ) );
%  layer indices of Green function objects
ind = arrayfun( @( x ) indlayer( x ), obj, 'uniform', 0 );
ind = vertcat( ind{ : } );
%  find Green function object for layer indices
[ ~, it ] = ismember( [ i1, i2 ], ind, 'rows' );
%  evaluate Green function object
switch p.Results.stat
  case 0
    y = eval( obj( it ), pos1, pos2, k0, varargin{ : } );
  case 1
    y = quasistatic( obj( it ), pos1, pos2, k0, varargin{ : } );
end
