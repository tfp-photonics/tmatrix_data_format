classdef tmatrix 
  %  T-matrix.
  
  properties
    aa      %  TM-TM elements
    ab      %  TM-TE elements
    ba      %  TE-TM elements
    bb      %  TE-TE elements
    k0      %  wavenumber of light in vacuum
  end
  
  properties (Hidden)
    solver  %  T-matrix or Mie solver
  end
  
  properties (Hidden, Dependent)
    k       %  wavenumber in embedding medium
    tab     %  table of spherical orders and degrees
  end  
  
  methods 
    function obj = tmatrix( varargin )
      %  Initialize T-matrix.
      %
      %  Usage :
      %    obj = multipole.tmatrix( solver, k0 )
      %  Input
      %    tsolver  :  T-matrix or Mie solver
      %    k0       :  wavenumber of light in vacuum
      if ~isempty( varargin ), obj = init( obj, varargin{ : } );  end
    end
    
    function sol = mldivide( obj, q )
      %  Solve T-matrix equations.
      sol = solve( obj, q );
    end       

    function x = get.k( obj )
      %  Wavenumber in embedding medium.
      x = obj.solver.embedding.k( obj.k0 );
    end    
    
    function x = get.tab( obj )
      %  Table of spherical orders and degrees.
      x = obj.solver.tab;
    end    
    
    function csca = scattering( obj )
      %  Average scattering cross section.
      csca = 2 * pi / obj.k ^ 2 * sum( abs( full( obj ) ) .^ 2, 'all' );
    end    
    
    function cext = extinction( obj )
      %  Average extinction cross section.
      cext = - 2 * pi / obj.k ^ 2 * real( trace( obj.aa ) + trace( obj.bb ) );
    end
    
   function cabs = absorption( obj )
      %  Average absorption cross section.
      cabs = extinction( obj ) - scattering( obj ); 
    end         
  end
end 
