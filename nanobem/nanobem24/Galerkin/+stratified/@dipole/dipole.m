classdef dipole
  %  Dipole excitation for stratified media.
  
  properties
    layer   %  layer structure
    pt      %  dipole positions
  end
  
  properties (Hidden)
    dip     %  galerkin.dipole object
    opt     %  options for reflected potential integrator
  end
    
  methods 
    function obj = dipole( varargin )
      %  Initialize dipole object.
      %
      %  Usage :
      %    obj = stratified.dipole( layer, pt, varargin )
      %  Input
      %    layer  :  layer structure
      %    pt     :  dipole positions
      obj = init( obj, varargin{ : } );
    end
  end

end 
