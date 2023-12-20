function varargout = eval( varargin )
%  EVAL - Evaluate Feibelman parameters.
%
%  Usage for obj = feibelman.param :
%    obj = eval( obj, k0 )
%    obj = eval( obj, tau, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum
%    tau    :  boundary elements

switch nargin
  case 2
    [ varargout{ 1 : nargout } ] = eval1( varargin{ : } );
  case 3
    [ varargout{ 1 : nargout } ] = eval2( varargin{ : } );
end
