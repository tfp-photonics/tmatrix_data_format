function obj = select( obj, fun )
%  SELECT - Select specific T-matrix elements.
%
%  Usage for obj = multipole.tmatrix :
%    t = select( obj, fun )
%  Input
%    fun  :  function handle @(m,l) or index to selected elements 
%  Output
%    obj  :  T-matrix with selected elements

%  deal with function handle
switch class( fun )
  case 'function_handle'
    tab = obj( 1 ).solver.tab;
    ind = fun( tab.m, tab.l );
  otherwise
    ind = fun;
end

switch numel( obj )
  case 1
    obj.aa = obj.aa( ind, ind );
    obj.ab = obj.ab( ind, ind );
    obj.ba = obj.ba( ind, ind );
    obj.bb = obj.bb( ind, ind );
  otherwise
    obj = arrayfun( @( x ) select( x, ind ), obj, 'uniform', 1 );
end
