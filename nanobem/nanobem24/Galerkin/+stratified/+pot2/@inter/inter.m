classdef inter < stratified.pot2.base
  %  Reflected interlayer BEM potentials.
  
  methods
    function obj = inter( varargin )
      %  Initialize reflected intralayer potential evaluator.
      %
      %  Usage :
      %    obj = stratified.pot2.inter( layer, tau1, tau2, PropertyPairs )
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
