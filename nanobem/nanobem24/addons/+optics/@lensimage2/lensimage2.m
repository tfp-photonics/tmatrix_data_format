classdef lensimage2
  %  Imaging lens using FFT, Hohenester Eq. (3.10).  
  
  properties
    mat1        %  material properties on object side
    mat2        %  material properties on image  side
    k0          %  wavenumber of light in vacuum
    NA          %  numerical aperture
    rot         %  rotation matrix for optical axis
    dir         %  directions for far-fields
    backfocal   %  field manipulation in backfocal plane
  end
  
  properties (Hidden)
    ind         %  index to directions inside of cutoff angle
    theta       %  cutoff angle for sphere at infinty
  end
  
  methods
    function obj = lensimage2( varargin )
      %  Initialize imaging lens using FFT.
      %
      %  Usage :
      %    obj = optics.lensimage2( mat1, mat2, k0, NA, PropertyPair )
      %  Input
      %    mat1       :  material properties on object side
      %    mat2       :  material properties on image  side
      %    k0         :  wavenumber of light in vacuum
      %    NA         :  numerical aperture
      %  PropertyName
      %    n          :  number of wavenumbers per dimension      
      %    rot        :  rotation matrix for optical axis
      %    backfocal  :  field manipulation in backfocal plane      
      obj = init( obj, varargin{ : } );
    end
  end
  
end
