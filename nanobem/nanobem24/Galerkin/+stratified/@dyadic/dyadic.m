classdef dyadic
  %  Dyadic Green function.
  
  properties
    tab       %  tabulated Green functions
  end
  
  methods
    function obj = dyadic( tab )
      %  Initialize dyadic Green function.
      %
      %  Usage :
      %    obj = stratified.dyadic( tab )
      %  Input
      %    tab    :  tabulated Green function object
      obj.tab = tab;
    end  
  end
end
