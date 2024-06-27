classdef bemsolver
  %  BEM solver for full Maxwell equations.
  
  properties
    tau     %  boundary elements
    pot     %  potential integrator
  end
  
  properties (Hidden = true)
    cal = []    %  Calderon matrix, for storage mode only
    k0          %  wavenumber for which Calderon matrix has been computed
  end
  
  methods 
    function obj = bemsolver( varargin )
      %  Initialize BEM solver for full Maxwell equations.
      %
      %  Usage :
      %    obj = galerkin.bemsolver( tau, PropertyPairs )
      %  Input
      %    tau        :  boundary elements
      %  PropertyName
      %    'engine'   :  potential integrator
      obj = init( obj, varargin{ : } );
    end
    
    function sol = mldivide( obj, q )
      %  MLDIVIDE - Solve BEM equations.
      sol = solve( obj, q );
    end 
    
    function cal = calderon( obj, k0 )
      %  CALDERON - Calderon matrix.
      cal = calderon( obj.pot, k0 );
    end    
  end
end  
