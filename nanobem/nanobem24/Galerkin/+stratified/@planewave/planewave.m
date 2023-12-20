classdef planewave
  %  Plane wave excitation for full Maxwell equations and layer structure.
  
  properties
    layer   %  layer structure
    pol     %  polarizations of plane waves
    dir     %  light propagation direction
  end
  
  properties (Hidden)
    rules   %  quadrature rules 
  end
  
  methods 
    function obj = planewave( varargin )
      %  Initialize plane wave object for layer structure.
      %
      %  Usage :
      %    obj = stratified.planewave( layer, pol, dir, varargin )
      %  Input
      %    layer  :  layer structure
      %    pol    :  polarizations of plane waves
      %    dir    :  propagation direction
      obj = init( obj, varargin{ : } );
    end
  end

end 
