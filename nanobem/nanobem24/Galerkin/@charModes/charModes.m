classdef charModes
  %  Characteristic modes for BEM solver.
  
  properties
    bem     %  BEM solver
    val     %  eigenvalues
    vec     %  eigenvectors
    k0      %  wavenumber of light in vacuum
    rsym    %  symmetry number, zero for symmetric Calderon matrix
  end
  
  properties (Hidden)
    A       %  symmetrized transmission matrix
  end
  
  
  methods
    function obj = charModes( bem )
      %  Initialize characteristic modes for BEM solver.
      %
      %  Usage :
      %    obj = charModes( bem )
      %  Input
      %    bem    :  BEM solver
      obj.bem = bem;
    end
  end
end
