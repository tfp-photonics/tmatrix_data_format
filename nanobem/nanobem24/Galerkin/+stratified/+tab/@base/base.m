classdef (Abstract) base  < matlab.mixin.Heterogeneous
  %  Base class for tabulated Green functions.
    
  properties
    layer     %  layer structure
    k0        %  wavenumber for tabulated Green functions
    method    %  interpolation method
  end
  
  properties (Hidden)
    somint    %  Sommerfeld integrator
  end

  methods
    function obj = base( varargin )
      %  Initialize base class for tabulated Green functions.
      %
      %  Usage :
      %    obj = stratified.tab.base( layer, i1, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    i1       :  medium index
      %  PropertyPairs for layerstructure.isommerfeld
      obj = init1( obj, varargin{ : } );
    end
  end
  
  methods (Sealed)
    function obj = fill( obj, k0 )
      %  FILL - Fill Green function tables.
      if isempty( obj( 1 ).k0 ) || obj( 1 ).k0 ~= k0
        obj = arrayfun( @( x ) fill2( x, k0 ), obj, 'uniform', 0 );
        obj = horzcat( obj{ : } );
      end
    end
    
    y = interp(  obj, varargin );
  end
end
