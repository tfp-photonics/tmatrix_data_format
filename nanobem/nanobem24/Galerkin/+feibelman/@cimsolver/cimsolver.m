classdef cimsolver < galerkin.cimsolver
  %  Solver for contour integral method using Feibelman parameters.
  
  properties
    param     %  Feibelman parameters
  end
     
  methods
    function obj = cimsolver( varargin )
      %  Initialize contour integral method solver.
      %
      %  Usage :
      %    obj = galerkin.cimsolver( tau, param, PropertyPairs )
      %  Input
      %    tau      :  boundary elements
      %    param    :  Feibelman parameters
      %  PropertyName
      %    contour  :  energies and weights for contour integration
      %    nr       :  number of columns for random matrix
      %    nz       :  maximal power for z^n
      %    seed     :  seed for random number generator
      obj = obj@galerkin.cimsolver( varargin{ 1 }, varargin{ 3 : end } );
      obj.param = varargin{ 2 };
    end
    
    function sol = match( obj, varargin )
      %  MATCH - Match tangential fields from boundary outside to inside.
      sol = match( feibelman.bemsolver(  ...
        obj.tau, obj.param, 'order', [] ), varargin{ : } );
    end
  end
end
