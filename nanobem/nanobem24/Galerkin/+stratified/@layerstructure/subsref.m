function varargout = subsref( obj, s )
%  Access to class functions and properties of layer structure.
    
switch s( 1 ).type  
  case '.'  
    switch s( 1 ).subs
      case 'kz'
        %  extract input
        [ k0, kpar, i1 ] = deal( s( 2 ).subs{ : } );
        %  evaluate material parameters
        obj = eval( obj, k0 );
        data = obj.data;
        %  z-component of wavevector
        varargout{ 1 } = stratified.zsqrt(  ...
                 data.eps( i1 ) * data.mu( i1 ) * k0 ^ 2 - kpar .^ 2 );
      otherwise
        [ varargout{ 1 : nargout } ] = builtin( 'subsref', obj, s );
    end
end
