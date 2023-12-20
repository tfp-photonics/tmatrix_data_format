classdef inter < stratified.pot1.base
  %  Reflected interlayer BEM potentials.
  
  methods
    function obj = inter( varargin )
      %  Initialize reflected interlayer potential evaluator.
      %
      %  Usage :
      %    obj = stratified.pot1.inter( layer, tau1, tau2, PropertyPairs )
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
