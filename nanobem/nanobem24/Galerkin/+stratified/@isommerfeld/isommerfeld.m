classdef isommerfeld
  %  Sommerfeld integrator.
  
  properties
    layer     %  layer structure
    i1        %  material index
    semi1     %  scaling factor for real axis of semiellipse
    semi2     %  scaling factor for imaginary axis of semiellipse
    ratio     %  z : r ratio which determines integration path
    op        %  options for ODE integration    
    cutoff    %  cutoff parameter for matrix-friendly approach of Chew
  end
  
  methods
    function obj = isommerfeld( varargin )
      %  Initialize Sommerfeld integrator.
      %
      %  Usage :
      %    obj = stratified.isommerfeld( layer, i1, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    i1       :  material index
      %  PropertyName
      %    semi1    :  scaling factor for real axis of semiellipse
      %    semi2    :  scaling factor for imaginary axis of semiellipse 
      %    ratio    :  z : r ratio which determines integration path
      %    op       :  options for ODE integration      
      %    cutoff   :  cutoff parameter for matrix-friendly approach of Chew
      obj = init( obj, varargin{ : } );
    end
  end
  
end