classdef solution < multipole.base
  %  Multipole solution.
  
  properties
    a     %  TM scattering coefficients
    b     %  TE scattering coefficients
    ai    %  TM incoming coefficients
    bi    %  TE incoming coefficients
    k0    %  wavenumber of light in vacuum
    mat   %  material properties of embedding medium
  end
  
  methods 
    function obj = solution( lmax, varargin )
      %  Initialize multipole solution.
      %
      %  Usage :
      %    obj = multipole.solution( lmax, mat, k0, a, b, ai, bi )
      %  Input
      %    lmax     :  maximal degree for multipole expansion
      %    mat      :  material properties of embedding medium
      %    k0       :  wavenumber of light in vacuum
      %    a        :  TM scattering coefficients
      %    b        :  TE scattering coefficients
      %    ai       :  TM incoming coefficients
      %    bi       :  TE incoming coefficients      
      obj = obj@multipole.base( lmax );
      obj = init( obj, varargin{ : } );
    end
    
    function cext = extinction( obj )
      %  Extinction cross section.
      fac = 0.5 * obj.mat.Z( obj.k0 ) / obj.mat.k( obj.k0 ) ^ 2;
      cext = - fac * real( obj.ai' * obj.a + obj.bi' * obj.b ); 
    end   
    
   function cabs = absorption( obj )
      %  Absorption cross section.
      cabs = extinction( obj ) - scattering( obj ); 
    end        
  end
end 
