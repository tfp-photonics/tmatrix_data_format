classdef bemsolver
  %  BEM solver for Maxwell equations and layer structure.
  
  properties
    tau     %  boundary elements
    layer   %  layer structure
    pot1    %  direct potential integrator
    pot2    %  reflected potential integrator
    green   %  tabulated Green function object
  end
  
  properties (Hidden = true)
    cal = []    %  LU decomposition of Calderon matrix
    k0          %  wavenumber for which Calderon matrix has been computed
  end
  
  methods 
    function obj = bemsolver( varargin )
      %  Initialize BEM solver for layer structures.
      %
      %  Usage :
      %    obj = stratified.bemsolver( tau, layer, PropertyPairs )
      %  Input
      %    tau        :  boundary elements
      %  PropertyName
      %    'engine'   :  engine for direct potential integrator
      obj = init( obj, varargin{ : } );
    end
    
    function sol = mldivide( obj, q )
      %  MLDIVIDE - Solve BEM equations.
      sol = solve( obj, q );
    end     
  end
end  
