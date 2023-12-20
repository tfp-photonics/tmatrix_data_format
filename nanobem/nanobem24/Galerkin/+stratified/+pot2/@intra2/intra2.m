classdef intra2 < stratified.pot2.base
  %  Reflected BEM potentials in given layer.
  
  methods
    function obj = intra2( varargin )
      %  Initialize reflected intralayer potential evaluator.
      %
      %  Usage :
      %    obj = stratified.pot2.intra2( layer, tau1, tau2, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    tau1     :  evaluation points
      %    tau2     :  boundary elements
      %  PropertyName
      %    rules    :  default quadrature rules
      obj = obj@stratified.pot2.base( varargin{ : } );
    end  
  end
end
