function varargout = subsref( obj, s )
%  Evaluate T-matrix solver.

switch s( 1 ).type
  case '()'
    varargout{ 1 } = qinc( obj, s( 1 ).subs{ : } );
  otherwise
    [ varargout{ 1 : nargout } ] = builtin( 'subsref', obj, s );
end
