classdef lensimage
  %  Imaging lens, Hohenester "Nano and Quantum Optics", Eq. (3.10).  
  
  properties
    mat1        %  material properties on object side
    mat2        %  material properties on image  side
    k0          %  wavenumber of light in vacuum
    NA          %  numerical aperture
    mcut        %  cutoff for angular degree
    rot         %  rotation matrix for optical axis
    dir         %  light propagation directions
    backfocal   %  field manipulation in backfocal plane
  end
  
  properties (Hidden)
    phi         %  azimuthal angles for sphere at infinity
    theta       %  polar angles for sphere at infinty
    w           %  weights for polar angle integration
  end
  
  methods
    function obj = lensimage( varargin )
      %  Initialize imaging lenses.
      %
      %  Usage :
      %    obj = optics.lensimage( mat1, mat2, k0, NA, PropertyPair )
      %  Input
      %    mat1       :  material properties on object side
      %    mat2       :  material properties on image  side
      %    k0         :  wavenumber of light in vacuum
      %    NA         :  numerical aperture
      %  PropertyName
      %    rot        :  rotation matrix for optical axis
      %    backfocal  :  field manipulation in backfocal plane
      %    nphi       :  number of azimuthal angles
      %    ntheta     :  number of polar angles
      %    mcut       :  cutoff for angular degree
      obj = init( obj, varargin{ : } );
    end
  end
  
end
