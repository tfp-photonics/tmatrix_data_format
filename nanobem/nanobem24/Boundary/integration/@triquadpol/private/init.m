function obj = init( obj, tau, pos, varargin )
%  INIT - Initialize polar triangle integration.
%
%  Usage for obj = triquadpol :
%    obj = init( obj, tau, pos, PropertyPairs )
%  Input
%    tau    :  boundary elements
%    pos    :  origin positions for polar integration
%  PropertyName
%    npol   :  number of radial and azimuthal integration points  

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'npol', [ 5, 5 ] );
%  parse input
parse( p, varargin{ : } );

%  assert triangular boundary elements 
n = unique( horzcat( tau.nedges ) );
assert( isscalar( n ) && n == 3 );
%  save input and adapt integration rules to boundaries
[ obj.tau, obj.pos, obj.npol ] = deal( tau, pos, p.Results.npol ); 
obj = adapt( obj );
