classdef base
  %  Table of spherical degrees and orders
  
  properties
    lmax    %  maximal degree for multipole expansion
    tab     %  table of spherical degrees and orders
  end
  
  methods 
    function obj = base( varargin )
      %  Initialize spherical degrees and orders.
      %
      %  Usage :
      %    obj = multipole.base( lmax )
      %    obj = multipole.base( tab )
      %  Input
      %    lmax     :  maximal degree for multipole expansion
      %    tab      :  struct with spherical degrees and orders
      obj = init( obj, varargin{ : } );
    end
  end

end 
