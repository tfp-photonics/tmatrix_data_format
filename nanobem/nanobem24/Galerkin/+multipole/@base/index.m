function ind = index( obj, varargin )
%  INDEX - Combined index for angular degrees and orders.
%
%  Usage for obj = multipole.base :
%    ind = index( obj, param )
%  Input
%    param  :  additional parameters for ind(:,3)

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'param', [] );
%  parse input
parse( p, varargin{ : } );

%  combined indes for angular degrees and orders
ind = [ obj.tab.l( : ), obj.tab.m( : ) ];
%  deal with additional parameter
if ~isempty( p.Results.param )
  fun = @( x ) [ ind, repelem( x, size( ind, 1 ), 1 ) ];
  ind = arrayfun( fun, p.Results.param, 'uniform', 0 );
  ind = vertcat( ind{ : } );
end
