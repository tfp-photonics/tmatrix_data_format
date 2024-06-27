function [ e, h ] = farfields( obj, dir, k0, varargin )
%  FARFIELDS - Electromagnetic farfields.
%
%  Usage for obj = laguerregauss :
%    [ e, h ] = farfields( obj, pos, k0, PropertyPairs )
%  Input
%    dir    :  propagation directions of farfields
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    infty  :  large number
%    shift  :  shift farfields
%  Output
%    e,h    :  electromagnetic farfields

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'infty', 1e10 );
addParameter( p, 'shift', [] );
%  parse input
parse( p, varargin{ : } );

%  paraxial fields at large distances
[ e, h ] = paraxial( obj, p.Results.infty * dir, k0 );
%  convert to farfield amplitudes
fac = exp( 1i * obj.mat.k( k0 ) * p.Results.infty ) / p.Results.infty;
[ e, h ] = deal( e / fac, h / fac );

%  deal with additional shift of positions
if ~isempty( p.Results.shift )
  %  phase factor
  fac = exp( - 1i * k * dir * p.Results.shift .' );
  %  multiply farfields with phase factor
  e = bsxfun( @times, e, fac );
  h = bsxfun( @times, h, fac );
end
