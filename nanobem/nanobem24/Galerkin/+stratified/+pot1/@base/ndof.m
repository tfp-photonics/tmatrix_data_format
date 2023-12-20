function n = ndof( obj, varargin )
%  NDOF - Global degrees of freedom.
%
%  Usage for obj = stratified.pot1.base :
%    n = ndof( obj, i1 )
%  Input
%    i1     :  global degrees of freedom for (1) tau1 or (2) tau2
%  Output
%    n      :  global degrees of freedom

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'i1', 1 );
parse( p, varargin{ : } );

switch p.Results.i1
  case 1
    n = max( arrayfun( @( x ) ndof( x.tau1 ), obj, 'uniform', 1 ) );
  case 2
    n = max( arrayfun( @( x ) ndof( x.tau2 ), obj, 'uniform', 1 ) );
end
