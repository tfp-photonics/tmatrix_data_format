classdef capillary
  %  Scattered light through capillary, Hohenester et al.,  13, 457 (2024).
    
  properties
    mat        %  material properties of stratified cylinders
    diameter   %  diameters of stratified cylinders
  end
    
  methods
    function obj = capillary( varargin )
      %  Initialize capillary.
      %
      %  Usage :
      %    obj = optics.capillary( mat, diameter )
      %  Input
      %    mat        :  material properties of stratified cylinders
      %    diameter   :  diameters of stratified cylinders
      if nargin,  [ obj.mat, obj.diameter ] = deal( varargin{ : } );  end
    end
  end  
end
