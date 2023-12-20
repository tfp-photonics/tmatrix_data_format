classdef solution < galerkin.solution
  %  Solution of BEM equation. for layer structure.
  
  properties
    layer   %  layer structure
  end
  
  methods
    function obj = solution( varargin )
      %  Initialize BEM solution for layer structure.
      %
      %  Usage :
      %    obj = solution( tau, k0, e, h )
      %  Input
      %    tau    :  vector of boundary elements
      %    k0     :  wavelength of light in vacuum
      %    e      :  electric field coefficients
      %    h      :  magnetic field coefficients
      obj = obj@galerkin.solution( varargin{ : } );
    end
  end 
end
