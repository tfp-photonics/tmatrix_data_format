classdef decompose
  %  Planewave decomposition.  
  
  properties
    k0      %  wavenumber of light in vacuum
    efield  %  electric field components
    dir     %  propagation directions
  end
  
  properties (Dependent, Hidden)
    n       %  number of field components
  end
  
  methods
    function obj = decompose( varargin )
      %  Initialize plane wave decomposition.
      %
      %  Usage :
      %    obj = decompose( k0, efield, dir )
      %  Input
      %    k0       :  wavenumber of light in vacuum
      %    efield   :  electric field components
      %    dir      :  propagation directions
      obj = init( obj, varargin{ : } );
    end
    
    function n = get.n( obj )
      n = size( obj.efield, 1 );
    end
  end
end
