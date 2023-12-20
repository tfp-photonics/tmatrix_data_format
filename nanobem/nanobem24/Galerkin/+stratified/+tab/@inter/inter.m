classdef inter < stratified.tab.base
  %  Tabulated Green functions connecting different layers.
    
  properties
    i1        %  medium index for observation points
    i2        %  medium index for source points
    rtab      %  radii for tabulated Green functions
    ztab1     %  z-values observation points of tabulated Green functions 
    ztab2     %  z-values source points of tabulated Green functions 
  end
  
  properties (Hidden)
    ytab      %  tabulated Green function
  end
  
  methods
    function obj = inter( layer, i1, i2, varargin )
      %  Initialize tabulated Green functions connecting different layers.
      %
      %  Usage :
      %    obj = stratified.tab.inter(  ...
      %            layer, i1, i2, rtab, ztab1, ztab2, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    i1       :  medium index for observation points
      %    i2       :  medium index for source points
      %    rtab     :  radii for tabulated Green functions
      %    ztab1    :  z-values observation points of tabulated Green functions 
      %    ztab2    :  z-values source points of tabulated Green functions 
      %  PropertyPairs for stratified.isommerfeld  
      obj = obj@stratified.tab.base( layer, i1, varargin{ 4 : end } );
      obj = init( obj, i1, i2, varargin{ : } );
    end
    
    function ind = indlayer( obj )
      %  INDLAYER - Index to layers.
      ind = [ obj.i1, obj.i2 ];
    end    
  end
end
