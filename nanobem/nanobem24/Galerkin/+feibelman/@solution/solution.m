classdef solution < galerkin.solution
  %  Solution of BEM equation with Feibelman parameters.
  
  properties 
    inout = 2     %  fields E,H at inside or outside
  end
  
  properties (Hidden)
    et      %  electric field at particle outside or inside
    ht      %  magnetic field at particle outside or inside
  end
  
  methods
    function obj = solution( varargin )
      %  Initialize BEM solution.
      %
      %  Usage :
      %    obj = feibelman.solution( tau, k0, e, h )
      %  Input
      %    tau    :  vector of boundary elements
      %    k0     :  wavelength of light in vacuum
      %    e      :  electric field coefficients
      %    h      :  magnetic field coefficients
      obj = obj@galerkin.solution( varargin{ : } ); 
    end
  end 
end
