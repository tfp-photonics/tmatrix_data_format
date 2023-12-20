classdef tsolver < multipole.base
  %  T-matrix solver.
 
  %  We use Eq. (E.4) of Hohenester "Nano and Quantum Optics" for the
  %  scattered and incoming fields.
  
  properties
    mat     %  material vector
    imat    %  index for embedding medium
  end
  
  properties (Hidden)
    rules       %  quadrature rules 
  end
  
  properties (Dependent, Hidden)
    embedding   %  material properties of embedding medium
  end  
  
  methods 
    function obj = tsolver( mat, imat, lmax, varargin )
      %  Initialize T-matrix solver.
      %
      %  Usage :
      %    obj = multipole.tsolver( mat, imat, lmax, PropertyPairs )
      %  Input
      %    mat      :  material vector
      %    imat     :  index for embedding medium      
      %    lmax     :  maximal degree for multipole expansion
      %  PropertyName
      %    rules    :  quadrature rules
      obj = obj@multipole.base( lmax );
      obj = init( obj, mat, imat, varargin{ : } );
    end
    
    function x = get.embedding( obj )
      %  EMBEDDING - Material properties of embedding medium.
      x = obj.mat( obj.imat );
    end
    
  end
end 
