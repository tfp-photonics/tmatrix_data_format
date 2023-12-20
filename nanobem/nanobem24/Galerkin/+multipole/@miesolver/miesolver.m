classdef miesolver < miesolver 
  %  Mie solver for multipole expansion, mainly implemented for testing.
    
  properties
    tab     %  table of spherical orders and degrees
  end
  
  properties (Dependent, Hidden)
    embedding   %  material properties of embedding medium
    mat         %  material vector
  end   
  
  properties (Hidden)
    imat = 1    %  index to embedding medium
  end     
  
  methods
    function obj = miesolver( mat1, mat2, diameter, varargin )
      %  Set up Mie solver.
      %
      %  Usage :
      %    obj = multipole.miesolver( mat1, mat2, diameter, PropertyPair )
      %  Input
      %    mat1       :  material properties at sphere inside
      %    mat2       :  material properties at sphere outside
      %    diameter   :  diameter of sphere in nm    
      %  PropertyName
      %    lmax       :  maximum number of spherical degrees
      obj = obj@miesolver( mat1, mat2, diameter ); 
      obj = init( obj, varargin{ : } );
    end
    
    function x = get.embedding( obj )
      %  EMBEDDING - Material properties of embedding medium.
      x = obj.mat2;
    end    
    
    function x = get.mat( obj )
      %  MAT - Material properties.
      x = [ obj.mat2, obj.mat1 ];
    end        
  end
  
end
