classdef triquadpol
  %  Polar triangle integration.
  
  properties
    tau   %  boundary elements
    pos   %  origin positions for polar integration
    quad  %  quadrature rules
    i1    %  quadrature rules for tau(i1) and pos(i1)
  end
  
  properties (Hidden)  
    npol  %  number of radial and azimuthal integration points
  end
  
  properties (Hidden, Dependent)
    npts  %  number of integration points
  end
  
  methods 
    function obj = triquadpol( varargin )
      %  Initialize polar triangle integration.
      %
      %  Usage :
      %    obj = triquadpol( tau, pos, PropertyPairs )
      %  Input
      %    tau    :  boundary elements
      %    pos    :  origin positions for polar integration
      %  PropertyName
      %    npol   :  number of radial and azimuthal integration points 
      if ~isempty( varargin ),  obj = init( obj, varargin{ : } );  end
    end
    
    function n = get.npts( obj )
      %  Number of integration points
      n = numel( obj.quad.w );
    end
  end   
end
