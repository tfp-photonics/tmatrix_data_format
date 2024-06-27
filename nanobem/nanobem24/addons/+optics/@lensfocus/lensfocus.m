classdef lensfocus
  %  Focus lens, Hohenester "Nano and Quantum Optics", Eq. (3.20).  
  properties
    mat     %  material properties on focus side
    k0      %  wavenumber of light in vacuum 
    NA      %  numerical aperture
    rot     %  rotation matrix for optical axis
  end
  
  properties (Dependent, Hidden)
    x,y     %  coordinates transverse to optical axis
  end
  
  properties (Hidden)
    rad     %  cutoff radius for Gaussian sphere    
    rho     %  radius, incoming field before Gaussian sphere
    phi     %  azimuthal angle   
    theta   %  polar angle, field after Gaussian sphere
    w       %  integration weights
  end
  
  methods
    function obj = lensfocus( varargin )
      %  Initialize imaging lenses.
      %
      %  Usage :
      %    obj = optics.lensfocus( mat, k0, NA, PropertyPair )
      %  Input
      %    mat      :  material properties on focus side
      %    k0       :  wavenumber of light in vacuum
      %    NA       :  numerical aperture
      %  PropertyName
      %    nphi     :  number of azimuthal angles
      %    ntheta   :  number of polar angles
      %    rot      :  rotation matrix for optical axis
      obj = init( obj, varargin{ : } );
    end
    
    function x = get.x( obj )
      x = cos( obj.phi ) .* obj.rho;
    end
    
    function y = get.y( obj )
      y = sin( obj.phi ) .* obj.rho;
    end    
    
  end
end
