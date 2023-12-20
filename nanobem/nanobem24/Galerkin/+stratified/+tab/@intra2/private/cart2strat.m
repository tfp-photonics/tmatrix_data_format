function [ r, Z1, Z2, R1, R2, siz ] = cart2strat( obj, pos1, pos2, varargin )
%  CART2STRAT - Convert to coordinates for stratified system.
%
%  Usage for obj = stratified.tab.intra2 :
%    [ r, Z1, Z2, R1, R2, siz ] = cart2strat( obj, pos1, pos2, PropertyPairs )
%  Input
%    pos1     :  observation points
%    pos2     :  source points
%  PropertyName
%    rows     :  positions of same size or not
%  Output
%    r        :  radial distances
%    Z1       :  distance from lower interface
%    Z2       :  distance from upper interface
%    R1       :  sqrt( r ^ 2 + Z1 ^ 2 )
%    R2       :  sqrt( r ^ 2 + Z2 ^ 2 )
%    siz      :  output size of arrays

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rows', 0 );
%  parse input
parse( p, varargin{ : } );

%  reshape position arrays
[ siz1, pos1 ] = deal( size( pos1 ), reshape( pos1, [], 3 ) );
[ siz2, pos2 ] = deal( size( pos2 ), reshape( pos2, [], 3 ) );
%  lower and upper interface position
z1 = obj.layer.z( obj.i1 - 1 );
z2 = obj.layer.z( obj.i1 );

switch p.Results.rows
  case 0
    %  radial distance and z-value
    r  = pdist2( pos1( :, 1 : 2 ), pos2( :, 1 : 2 ) );
    Z1 = abs( bsxfun( @plus, pos1( :, 3 ), pos2( :, 3 ) .' ) - 2 * z1 );
    Z2 = abs( bsxfun( @plus, pos1( :, 3 ), pos2( :, 3 ) .' ) - 2 * z2 );
    %  size of output array
    siz = [ siz1( 1 : end - 1 ), siz2( 1 : end - 1 ) ];
  case 1
    %  radial distance and z-value
    r  = sqrt( ( pos1( :, 1 ) - pos2( :, 1 ) ) .^ 2 +  ...
               ( pos1( :, 2 ) - pos2( :, 2 ) ) .^ 2 );
    Z1 = abs( pos1( :, 3 ) + pos2( :, 3 ) - 2 * z1 );
    Z2 = abs( pos1( :, 3 ) + pos2( :, 3 ) - 2 * z2 );
    %  size of output array
    siz = siz1( 1 : end - 1 );
end

%  full distance
R1 = sqrt( r .^ 2 + Z1 .^ 2 );
R2 = sqrt( r .^ 2 + Z2 .^ 2 );
