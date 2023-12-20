function y = green( grid, varargin )
%  GREEN - Initialize Green function objects.
%
%  Usage :
%    y = stratified.tab.green( grid, PropertyPairs )
%  Input
%    grid   :  tabulation grids, see stratified.layerstructure/grid
%  PropertyPairs to be passed to Sommerfeld integrator
%  Output
%    y      :  Green function tabulators

%  allocate output
y = cell( size( grid ) );
%  loop over tabulation grids
for it = 1 : numel( grid )
  x = grid( it );
  %  assign Green function objects
  if x.i1 == x.i2
    switch x.i1
      case { 1, x.layer.n + 1 }
        %  object for upper or lower layer
        y{ it } = stratified.tab.intra1(  ...
               x.layer, x.i1, x.r, x.z1, varargin{ : } );
      otherwise
        %  intralayer object
        y{ it } = stratified.tab.intra2(  ...
               x.layer, x.i1, x.r, x.z1, x.z2, varargin{ : } );
    end
  else
    %  interlayer object
    y{ it } = stratified.tab.inter(  ...
           x.layer, x.i1, x.i2, x.r, x.z1, x.z2, varargin{ : } );
  end
end

%  allocate output
y = horzcat( y{ : } );
