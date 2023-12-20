function varargout = subsref( obj, s )
%  Access to class functions and properties of Sommerfeld integrator.
    
switch s( 1 ).type  
  case '()'
    [ varargout{ 1 : nargout } ] = eval( obj, s( 1 ).subs{ : } );
  case '.'  
    [ varargout{ 1 : nargout } ] = builtin( 'subsref', obj, s );
end
