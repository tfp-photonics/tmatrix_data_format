function data = eval1( obj, varargin )
%  EVAL1 - Evaluate contour integrals.
%
%  Usage for obj = stratified.cimsolver :
%    data = eval1( obj,      PropertyPairs )
%    data = eval1( obj, bem, PropertyPairs )
%  Input
%    bem    :  boundary element method solver

if ~isempty( varargin ) && ~ischar( varargin{ 1 } )
  data = eval1@cimbase( obj, varargin{ : } );
else
  bem = feibelman.bemsolver( obj.tau, obj.param, varargin{ : } );
  data = eval1@cimbase( obj, bem, varargin{ : } );
end
