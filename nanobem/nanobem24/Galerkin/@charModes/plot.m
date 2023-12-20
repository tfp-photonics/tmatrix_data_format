function plot( obj, varargin )
%  PLOT - Plot surface charge of characteristic modes.
%
%  Usage for obj = charModes :
%    plot( obj, PropertyPairs )
%  PropertyName
%    ind    :  index to selected modes

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'ind', 1 : 2 * ndof( obj.bem.tau ) );
%  parse input
parse( p, varargin{ : } );

assert( ~isempty( obj.vec ) );
%  characteristic modes
vec = obj.vec;
vec.e = vec.e( :, p.Results.ind );
vec.h = vec.h( :, p.Results.ind );
%  final plot
figure
plot( vec.tau, surfc( vec ) );
