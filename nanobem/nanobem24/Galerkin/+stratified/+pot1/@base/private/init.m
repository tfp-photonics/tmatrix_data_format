function obj = init( obj, layer, tau1, tau2, varargin )
%  INIT - Initialize reflected potential objects.
%
%  Usage for obj = stratified.pot.base :
%    obj = init( obj, layer, tau1, tau2, PropertyPairs )
%  Input
%    layer    :  layer structure
%    tau1     :  first set of boundary elements
%    tau2     :  second set of boundary elements
%  PropertyName
%    rules        :  default quadrature rules
%    rules1       :  refined quadrature rules
%    nduffy       :  number of Legendre-Gauss points for Duffy integration
%    layercutoff  :  relative cutoff for refined integration
%    memax        :  slice integration into bunches of size MEMAX 

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules',  quadboundary.rules( 'quad3', triquad(  3 ) ) );
addParameter( p, 'rules1', quadboundary.rules( 'quad3', triquad( 11 ) ) );
addParameter( p, 'nduffy', 3 );
addParameter( p, 'layercutoff', 2 );
addParameter( p, 'memax', 1e7 );
%  parse input
parse( p, varargin{ : } );

%  quadrature points
pt1 = quadboundary( tau1, p.Results.rules );  n1 = numel( pt1 );
pt2 = quadboundary( tau2, p.Results.rules );  n2 = numel( pt2 );

%  deal with multiple and single integration points
if n1 ~= 1 || n2 ~= 1
  %  allocate output
  obj = num2cell( repelem( obj, n1, n2 ) );
  %  allocate objects
  for i1 = 1 : n1
  for i2 = 1 : n2
    tau1 = horzcat( pt1( i1 ).tau );
    tau2 = horzcat( pt2( i2 ).tau );
    obj{ i1, i2 } = init( obj{ i1, i2 }, layer, tau1, tau2, varargin{ : } );
  end
  end
  %  assemble output
  obj = horzcat( obj{ : } );
  
else
  %  save layer structure and quadrature points
  [ obj.layer, obj.pt1, obj.pt2 ] = deal( layer, pt1, pt2 );
  %  touching elements
  [ d, iz ] = bdist2( obj.layer, tau1, tau2 );
  ind1 = touching( tau1, tau2 );
  %  allocate output
  obj.yout = [];

  for iz1 = unique( iz ) .'
    %  elements for refinement
    ind2 = d <= p.Results.layercutoff & iz == iz1;
    %  boundary element pairs with Duffy integration
    if nnz( ind1 & ind2 )
      %  boundary element pairs
      [ i1, i2 ] = find( ind1 & ind2 );
      %  Duffy integration points
      pts = quadduffy( tau1( i1 ), tau2( i2 ), p.Results.nduffy, 'rows', 1 );
      pts = slice1( pts, p.Results.memax );
      %  set output
      obj.yout = horzcat( obj.yout, struct( 'pts', num2cell( pts ), 'iz', iz1 ) );
    end
    %  boundary element pairs with refined integration
    if nnz( ~ind1 & ind2 )
      %  boundary element pairs
      [ i1, i2 ] = find( ~ind1 & ind2 ); 
      %  refined integration points
      pt1 = quadboundary( tau1( i1 ), p.Results.rules1 );
      pt2 = quadboundary( tau2( i2 ), p.Results.rules1 );
      pts = slice2( [ pt1, pt2 ], p.Results.memax );
      %  set output
      pts = mat2cell( pts, ones( 1, size( pts, 1 ) ), 2 );
      obj.yout = horzcat( obj.yout, struct( 'pts', pts .', 'iz', iz1 ) );    
    end
  end
end


function pts = slice1( pts, memax )
%  SLICE1 - Slice Duffy integration points.

switch numel( pts )
  case 1
    if stratified.pot1.npts( pts ) > memax
      pts = slice( pts );
      pts = slice1( pts, memax );
    end
  otherwise
    pts = arrayfun( @( x ) slice1( x, memax ), pts, 'uniform', 0 );
    pts = horzcat( pts{ : } );
end


function pts = slice2( pts, memax )
%  SLICE2 - Slice quadrature points into bunches of size MEMAX.

switch size( pts, 1 )
  case 1
    if stratified.pot1.npts( pts ) > memax
      pt1 = slice( pts( 1 ) );
      pt2 = slice( pts( 2 ) );
      pts = slice2(  ...
        [ pt1( [ 1, 1, 2, 2 ] ) .', pt2( [ 1, 2, 1, 2 ] ) .' ], memax );
    end
  otherwise
    pts = mat2cell( pts, ones( 1, size( pts, 1 ) ), 2 );
    pts = cellfun( @( x ) slice2( x, memax ), pts, 'uniform', 0 );
    pts = vertcat( pts{ : } );
end

