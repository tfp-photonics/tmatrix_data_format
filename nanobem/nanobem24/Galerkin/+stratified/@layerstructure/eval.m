function obj = eval( obj, k0 )
%  EVAL - Evaluate permittivities and permeabilities for different media.
%
%  Usage for obj = stratified.layerstructure :
%    obj = eval( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum

if isempty( obj.data ) || obj.data.k0 ~= k0
  %  material properties
  %    the evaluation of tabulated material properties may slow down the
  %    calculation of Sommerfeld integrals, so it is better to compute the
  %    permittivities and permeabilities only once
  [ eps, mu ] = arrayfun(  ...
    @( x ) deal( x.eps( k0 ), x.mu( k0 ) ), obj.mat, 'uniform', 1 );
  %  set material properties 
  obj.data = struct( 'k0', k0, 'eps', eps, 'mu', mu );
end
