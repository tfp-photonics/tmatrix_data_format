classdef intra1 < stratified.pot2.base
  %  Reflected BEM potentials for upper or lower medium.
  
  methods
    function obj = intra1( varargin )
      %  Initialize reflected intralayer potential evaluator.
      %
      %  Usage :
      %    obj = stratified.pot2.intra1( layer, tau1, tau2, PropertyPairs )
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
