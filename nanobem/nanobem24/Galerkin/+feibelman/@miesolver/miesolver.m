classdef miesolver < miesolver
  %  Mie theory using Feibelman parameters, Nat. Comm. 11, 366 (2020).
  
  properties
    dperp     %  perpendicular Feibelman parameter
    dpar      %  parallel Feibelman parameter
  end
  
  methods
    function obj = miesolver( varargin )
      %  Set up Mie solver.
      %
      %  Usage :
      %    obj = feibelman.miesolver( mat1, mat2, diameter, PropertyPair )
      %  Input
      %    mat1       :  material properties at sphere inside
      %    mat2       :  material properties at sphere outside
      %    diameter   :  diameter of sphere in nm    
      %  PropertyName
      %    lmax       :  maximum number of spherical degrees
      obj = obj@miesolver( varargin{ : } );    
    end
  end
  
end
