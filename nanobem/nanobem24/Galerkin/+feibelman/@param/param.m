classdef param
  %  Feibelman parameters.
  
  properties
    tau     %  boundary elements
    fun     %  evaluation function for Feibelman parameters
  end
  
  methods
    function obj = param( varargin )
      %  Set Feibelman parameters.
      %
      %  Usage :
      %    obj = feibelman.param( tau, fun )
      %  Input
      %    tau  :  boundary elements
      %    fun  :  evaluation function for Feibelman parameters
      [ obj.tau, obj.fun ] = deal( varargin{ : } );
    end
  end
end 
