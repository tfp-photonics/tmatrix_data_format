classdef bemsolver
  %  BEM solver for full Maxwell equations with Feibelman parameters.
  
  properties
    tau       %  boundary elements
    pot       %  potential integrator
    param     %  Feibelman parameters
    ind       %  connected edges
  end

  methods
    function obj = bemsolver( varargin )
      %  Initialize BEM solver with Feibelman parameters.
      %
      %  Usage :
      %    obj = feibelman.bemsolver( tau, param, PropertyPairs )
      %  Input
      %    tau    :  boundary elements
      %    param  :  Feibelman parameters
      obj = init( obj, varargin{ : } );
    end
    
    function varargout = mldivide( obj, q )
      %  MLDIVIDE - Solve BEM equations.
      [ varargout{ 1 : nargout } ] = solve( obj, q );
    end     
  end
end 
