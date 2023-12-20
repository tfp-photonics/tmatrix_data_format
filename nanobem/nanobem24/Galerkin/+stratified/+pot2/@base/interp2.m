function [ y1, y2 ] = interp2( ~, green, pos1, pos2, varargin )
%  INTERP1 - Interpolate quasistatic Green function and derivatives.
%
%  Usage for obj = stratified.pot2.base :
%    [ y1, y2 ] = interp2( obj, green, pos1, pos2, PropertyPairs )
%  Input
%    green  :  Green function object
%    pos1   :  observation points
%    pos2   :  source points
%  PropertyName
%    ind    :  tensor indices [i,q,k]
%  Output
%    y1     :  reflected Green functions
%    y2     :  gradient of reflected Green function

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ind', 1 : 3 );
%  parse input
parse( p, varargin{ : } );

eta = 1e-4;
%  expand observation points
pos1 = cat( 3, bsxfun( @plus, pos1, eta * [ 1, 0, 0 ] ),  ...
               bsxfun( @plus, pos1, eta * [ 0, 1, 0 ] ),  ...
               bsxfun( @plus, pos1, eta * [ 0, 0, 1 ] ), pos1 );
pos1 = permute( pos1, [ 1, 3, 2 ] );
%  bring POS1 and POS2 to same shape
siz = size( pos2 );
pos1 = repmat( reshape( pos1, siz( 1 ), 4, 1, 3 ), 1, 1, siz( 2 ), 1 );
pos2 = repmat( reshape( pos2, siz( 1 ), 1, siz( 2 ), 3 ), 1, 4, 1, 1 );
%  evaluate quasistatic Green functions
y = quasistatic( green, reshape( pos1, siz( 1 ) * 4, [], 3 ),  ...
                        reshape( pos2, siz( 1 ) * 4, [], 3 ), [], 'rows', 1 );

%  dummy indices for tensor class
ind = num2cell( p.Results.ind );
[ i, q, k ] = deal( ind{ : } );
%  allocate output
[ y1, y2 ] = deal( repelem( struct, 1, numel( y ) ) );

%  loop over reflected Green functions
for it = 1 : numel( y )
for name = convertCharsToStrings( fieldnames( y( it ) ) ) .'
  %  reflected Green function
  yy = reshape( y( it ).( name ), siz( 1 ), 4, siz( 2 ) );
  y1( it ).( name ) = reshape( yy( :, 4, : ), size( yy, 1 ), [] );
  y1( it ).( name ) = tensor( y1( it ).( name ), [ i, q ] );
  %  derivative of reflected Green function
  y2( it ).( name ) = cat( 2, yy( :, 1, : ) - yy( :, 4, : ),  ...
                              yy( :, 2, : ) - yy( :, 4, : ),  ...
                              yy( :, 3, : ) - yy( :, 4, : ) ) / eta;
  y2( it ).( name ) = tensor( y2( it ).( name ), [ i, k, q ] );                      
end
end
