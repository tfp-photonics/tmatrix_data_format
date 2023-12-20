function in = inside( obj, pos, varargin )
%  INSIDE - Locate positions inside of discretized particle boundary.
%
%  Usage for obj = BoundaryElement :
%    in = inside( obj, pos, PropertyPairs )
%  Input
%    pos    :  point positions
%  PropertyName
%    i1     :  material indices of embedding medium
%  Output
%    in     :  index to material within which POS is located

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'i1', 1 );
%  parse input
parse( p, varargin{ : } );

%  set material index to embedding medium
in = repelem( p.Results.i1( 1 ), size( pos, 1 ), 1 );
%  material tree
inout = vertcat( obj.inout );
t = tree( unique( inout, 'rows' ), p.Results.i1 );

%  loop over media
for imat = cat( 2, t{ 2 : end } )
  %  boundary elements connected to medium IMAT
  ind = any( bsxfun( @eq, vertcat( obj.inout ), imat ), 2 );
  %  unique vertices
  verts = round( vertices( obj( ind ) ), 5 );
  [ verts, ~, faces ] = unique( reshape( verts, [], 3 ), 'rows' );
  %  faces of boundary elements
  ind = inout( ind, 1 ) == imat;
  faces = reshape( faces, [], 3 );
  faces( ind, : ) = fliplr( faces( ind, : ) );
  %  find positions inside of boundary
  [ in( inpolyhedron( faces, verts, pos ) ) ] = deal( imat );
end


function t = tree( inout, i1 )
%  Tree - Material tree.

%  root of tree and unique materials
t{ 1 } = i1;
imat = setdiff( unique( inout ), i1 );

while ~isempty( imat )
  %  materials connected to outer materials
  ind = any( bsxfun( @eq, inout( :, 1 ), t{ end } ), 2 ) |  ...
        any( bsxfun( @eq, inout( :, 2 ), t{ end } ), 2 ) ;
  %  update tree
  t{ end + 1 } =  ...
    setdiff( reshape( unique( inout( ind, : ) ), 1, [] ), t{ end } );
  %  update auxiliary arrays
  inout = inout( ~ind, : );
  imat = setdiff( imat, t{ end } ); 
end
