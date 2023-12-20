classdef intra1 < stratified.pot1.base
  %  Reflected BEM potentials for upper or lower medium.
  
  methods
    function obj = intra1( varargin )
      %  Initialize reflected intralayer potential evaluator.
      %
      %  Usage :
      %    obj = stratified.pot1.intra1( layer, tau1, tau2, PropertyPairs )
      %  Input
      %    layer    :  layer structure
      %    tau1     :  first set of boundary elements
      %    tau2     :  second set of boundary elements
      %  PropertyName
      %    rules    :  default quadrature rules
      obj = obj@stratified.pot1.base( varargin{ : } );
    end  
  end
end
