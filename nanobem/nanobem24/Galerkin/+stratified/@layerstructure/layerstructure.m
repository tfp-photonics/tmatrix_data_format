classdef layerstructure
  %  Layer structure for stratified media.
  
  properties
    mat   %  table of material properties
    z     %  z-values of interfaces
  end
  
  properties (Hidden)
    d     %  layer thickness
    data  %  material data
  end
  
  properties (Dependent, Hidden)
    n     %  number of interfaces
  end
  
  methods
    function obj = layerstructure( varargin )
      %  Initialize layer structure.
      %
      %  Usage :
      %    obj = stratified.layerstructure( mat, z )
      %  Input
      %    mat  :  material properties
      %    z    :  z-values of interfaces
      obj = init( obj, varargin{ : } );
    end    
    
    function n = get.n( obj )
      %  Number of interfaces.
      n = numel( obj.z );
    end   
  end
end
