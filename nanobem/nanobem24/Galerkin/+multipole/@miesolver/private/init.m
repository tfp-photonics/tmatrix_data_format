function obj = init( obj, x )
%  INIT - Initialize Mie solver.
%
%  Usage :
%    obj = multipole.miesolver( lmax )
%    obj = multipole.miesolver( tab )
%  Input
%    lmax     :  maximal degree for multipole expansion
%    tab      :  struct with spherical degrees and orders

%  default value for maximal degree
if ~exist( 'x', 'var' ),  x = 5;  end

switch isnumeric( x )
  case 1
    obj.lmax = x;
    [ obj.tab.l, obj.tab.m ] = sphtable( x );
  otherwise
    obj.lmax = max( x.l );
    obj.tab = x;
end
