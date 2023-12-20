classdef plot
  %  Plot tabulated Green function.
  
  properties
    rtab      %  radii for tabulated Green functions
    ztab      %  z-values for tabulated Green functions    
    ytab      %  tabulated Green functions 
  end
  
  properties (Hidden)
    fun       %  @real, @imag, @abs, or user-defined function
    name      %  name of function to be plotted
  end
  
  methods
    function obj = plot( varargin )
      %  Initialize Green function plotter.
      %
      %  Usage :
      %    obj = stratified.plot( rtab, ztab, ytab, PropertyPair )
      %  Input
      %    rtab     :  radii for tabulated Green functions
      %    ztab     :  z-values for tabulated Green functions    
      %    ytab     :  tabulated Green function  
      %  PropertyPair
      %    fun      :  @real, @imag, @abs, or user-defined function
      obj = init( obj, varargin{ : } );
    end
  end
  
end