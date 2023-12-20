classdef (Abstract) base < matlab.mixin.Heterogeneous 
  %  Base class for evaluation of reflected BEM potentials.
    
  properties 
    layer     %  layer structure
    pt1       %  quadrature points for first boundary
    pt2       %  quadrature points for second boundary   
  end
  
  properties (Hidden)
    yout      %  quasistatic contributions with refinement
  end
  
  properties (Dependent, Hidden)
    tau1      %  first set of boundary elements
    tau2      %  second set of boundary elements    
    i1        %  material index for TAU1
    i2        %  material index for TAU2     
  end   
  
  methods
    function obj = base( varargin )
      %  Initialize reflected potential object.
      %
      %  Usage :
      %    obj = stratified.pot1.base( layer, tau1, tau2, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    tau1     :  first set of boundary elements
      %    tau2     :  second set of boundary elements
      %  PropertyName
      %    rules    :  default quadrature rules
      obj = init( obj, varargin{ : } );
    end
    
    %  boundary elements
    function tau = get.tau1( obj ),  tau = obj.pt1.tau;  end
    function tau = get.tau2( obj ),  tau = obj.pt2.tau;  end
    %  material indices
    function tau = get.i1( obj ),  tau = obj.pt1.inout( 2 );  end
    function tau = get.i2( obj ),  tau = obj.pt2.inout( 2 );  end     
  end  

  methods (Static)
    n = npts( varargin );
    data = waitbar( varargin );
  end
  
  %  seal methods for matlab.mixin.Heterogeneous objects
  methods (Sealed)
    varargout = calderon( varargin );
    varargout = eval( varargin );
    varargout = ndof( varargin );
    varargout = oint1( varargin );
    varargout = oint2( varargin );
  end  
end
