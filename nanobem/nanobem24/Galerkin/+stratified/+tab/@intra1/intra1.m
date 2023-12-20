classdef intra1 < stratified.tab.base
  %  Tabulated Green functions in upper or lower layer.
    
  properties
    i1        %  medium index
    rtab      %  radii for tabulated Green functions
    ztab      %  z-values for tabulated Green functions        
  end
  
  properties (Hidden)
    ytab      %  tabulated Green function
  end
  
  methods
    function obj = intra1( layer, i1, varargin )
      %  Initialize tabulated Green functions for upper or lower layer.
      %
      %  Usage :
      %    obj = stratified.tab.intra1( layer, i1, rtab, ztab, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    i1       :  layer index
      %    rtab     :  radii for tabulated Green functions
      %    ztab     :  z-values for tabulated Green functions    
      %  PropertyPairs for stratified.isommerfeld  
      obj = obj@stratified.tab.base( layer, i1, varargin{ 3 : end } );
      obj = init( obj, i1, varargin{ : } );
    end
    
    function ind = indlayer( obj )
      %  INDLAYER - Index to layers.
      ind = [ obj.i1, obj.i1 ];
    end
  end
end
