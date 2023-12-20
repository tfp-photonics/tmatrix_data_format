classdef intra2 < stratified.tab.base
  %  Tabulated Green functions within given layer.
    
  properties
    i1        %  medium index
    rtab      %  radii for tabulated Green function
    ztab1     %  z-values for tabulated Green function, z1-z2      
    ztab2     %  z-values for tabulated Green function, z1+z2
  end
  
  properties (Hidden)
    ytab1     %  tabulated Green function, z1-z2
    ytab2     %  tabulated Green function, z1+z2
  end
  
  methods
    function obj = intra2( layer, i1, varargin )
      %  Initialize tabulated Green functions ithin given layer.
      %
      %  Usage :
      %    obj = stratified.tab.intra2(  ...
      %                layer, i1, rtab, ztab1, ztab2, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    i2       :  layer index
      %    rtab     :  radii for tabulated Green functions
      %    ztab1    :  z-values for tabulated Green function, z1-z2
      %    ztab2    :  z-values for tabulated Green function, z1-z2
      %  PropertyPairs for stratified.isommerfeld  
      obj = obj@stratified.tab.base( layer, i1, varargin{ 4 : end } );
      obj = init( obj, i1, varargin{ : } );
    end
    
    function ind = indlayer( obj )
      %  INDLAYER - Index to layers.
      ind = [ obj.i1, obj.i1 ];
    end    
  end
end
