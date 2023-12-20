function obj = init( obj, layer, tau1, tau2, varargin )
%  INIT - Initialize reflected potential objects.
%
%  Usage for obj = stratified.pot.base :
%    obj = init( obj, layer, tau1, tau2, PropertyPairs )
%  Input
%    layer    :  layer structure
%    tau1     :  evaluation points
%    tau2     :  boundary elements
%  PropertyName
%    rules    :  default quadrature rules
%    npol     :  number of quadratue points for polar integration
%    memax    :  slice integration into bunches of size MEMAX 

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules',  quadboundary.rules( 'quad3', triquad( 3 ) ) );
addParameter( p, 'layercutoff', 2 );
addParameter( p, 'memax', 1e7 );
%  parse input
parse( p, varargin{ : } );

%  save layer structure
obj.layer = layer;
%  quadrature points, assert single quadrature object
obj.pt1 = tau1;
obj.pt2 = quadboundary( tau2, p.Results.rules );
assert( numel( obj.pt2 ) == 1 );

%  scaled distance between boundary elements and evaluation points 
[ d, iz ] = bdist2( obj.layer, tau2, vertcat( tau1.pos ) );
[ d, iz ] = deal( transpose( d ), transpose( iz ) );
%  allocate output
obj.yout = [];

for iz1 = unique( iz ) .'
  %  elements for refinement
  ind = d <= p.Results.layercutoff & iz == iz1;
  %  evaluation points and boundary elements with refined integration
  if nnz( ind ) > 1
    %  evaluation points and boundary elements
    [ i1, i2 ] = find( ind ); 
    %  refined polar integration rules
    pts = triquadpol( tau2( i2 ), tau1.pos( i1, : ), varargin{ : } );
    obj.yout = horzcat( obj.yout,  ...
      struct( 'pts', pts, 'i1', i1, 'i2', i2, 'iz', iz1 ) );  
  end
end

%  precompute quasistatic elements
if ~isempty( obj.yout ),  obj = refine1( obj );  end
  
