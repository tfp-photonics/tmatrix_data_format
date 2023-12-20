function [ e, h ] = interp( obj, varargin )
%  INTERP - Interpolate BEM solution from edges to faces.
%
%  Usage for obj = feibelman.solution :
%    [ e, h ] = interp( obj, pts, PropertyPairs )
%  Input
%    pts    :  quadrature points
%  PropertyName
%    inout  :  compute fields at boundary inside or outside
%  Output
%    e      :  electric fields at requested points
%    h      :  magnetic fields at requested points

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'pt', [] );
addParameter( p, 'inout', 2 );
%  parse input
parse( p, varargin{ : } );

[ e, h ] = interp@galerkin.solution( toggle( obj, p.Results.inout ),  ...
                          varargin{ : }, 'inout', p.Results.inout );
