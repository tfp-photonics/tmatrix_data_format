classdef intra2 < stratified.pot1.base
  %  Reflected BEM potentials within given layer.
  
  methods
    function obj = intra2( varargin )
      %  Initialize reflected intralayer potential evaluator.
      %
      %  Usage :
      %    obj = stratified.pot1.intra2( layer, tau1, tau2, PropertyPairs )
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
