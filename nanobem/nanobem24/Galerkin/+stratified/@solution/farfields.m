function [ e, h ] = farfields( obj, dir, varargin )
%  FARFIELDS - Electromagnetic farfields along requested directions.
%
%  Usage for obj = stratified.solution :
%    [ e, h ] = farfields( obj, dir, PropertyPairs )
%  Input
%    dir    :  light propagation directions
%  PropertyName
%    rules  :  quadrature rules
%  Output
%    e      :  electric farfields along requested directions
%    h      :  magnetic farfields along requested directions

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules', quadboundary.rules );
%  parse input
parse( p, varargin{ : } );

%  material parameters
data.eps = arrayfun( @( x ) x.eps( obj.k0 ), obj.layer.mat, 'uniform', 1 );
data.mu  = arrayfun( @( x ) x.mu ( obj.k0 ), obj.layer.mat, 'uniform', 1 );
data.Z   = arrayfun( @( x ) x.Z  ( obj.k0 ), obj.layer.mat, 'uniform', 1 );
data.k   = arrayfun( @( x ) x.k  ( obj.k0 ), obj.layer.mat, 'uniform', 1 );
%  boundary elements connected to layer structure
inout = vertcat( obj.tau.inout );
tau = obj.tau( inout( :, 2 ) <= obj.layer.n + 1 );

%  allocate output
siz = size( obj.e );
[ e, h ] = deal( zeros( size( dir, 1 ), 3, prod( siz ) / siz( 1 ) ) );
%  avoid directions parallel to z-direction
i1 = abs( dir( :, 3 ) ) == 1;
t = 1e-5;
dir( i1, : ) = dir( i1, : ) *  ...
  [ cos( t ), 0, sin( t ); 0, 1, 0; - sin( t ), 0, cos( t ) ];
%  downgoing and upgoing waves
i1 = dir( :, 3 ) < 0;
i2 = dir( :, 3 ) > 0;

%  loop over unique boundary elements 
for pt = quadboundary( tau, 'rules', p.Results.rules )  
  %  quadrature points, weights, and shape elements
  [ data.pos, data.w, data.f ] = eval( pt );
  %  tangential electromagnetic fields at particle boundary
  nu = vertcat( pt.tau.nu );
  ue = reshape( obj.e, size( obj.e, 1 ), [] );  data.ue = ue( nu( : ), : );
  uh = reshape( obj.h, size( obj.h, 1 ), [] );  data.uh = uh( nu( : ), : );
  %  reflected farfields, downgoing waves
  if nnz( i1 )
    [ e1, h1 ] = fields2( obj, data, dir( i1, : ), pt );
    e( i1, :, : ) = e( i1, :, : ) + e1;
    h( i1, :, : ) = h( i1, :, : ) + h1;
  end
  %  reflected farfields, upgoing waves
  if nnz( i2 )
    [ e2, h2 ] = fields2( obj, data, dir( i2, : ), pt );
    e( i2, :, : ) = e( i2, :, : ) + e2;
    h( i2, :, : ) = h( i2, :, : ) + h2;
  end 
end

%  direct contributions, downgoing waves
if nnz( i1 )
  pt1 = Point( obj.layer.mat, 1, dir( i1, : ) );
  [ e1, h1 ] = farfields@galerkin.solution( obj, pt1 );
  e( i1, :, : ) = e( i1, :, : ) + e1;
  h( i1, :, : ) = h( i1, :, : ) + h1;
end
%  direct contributions, upgoing waves  
if nnz( i2 )
  pt1 = Point( obj.layer.mat, obj.layer.n + 1, dir( i2, : ) );
  [ e2, h2 ] = farfields@galerkin.solution( obj, pt1 );
  e( i2, :, : ) = e( i2, :, : ) + e2;
  h( i2, :, : ) = h( i2, :, : ) + h2;
end

%  reshape output
siz = [ size( dir, 1 ), 3, siz( 2 : end ) ];
[ e, h ] = deal( reshape( e, siz ), reshape( h, siz ) );
