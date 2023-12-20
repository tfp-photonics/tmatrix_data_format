function y = eval( obj, pos1, pos2, k0, varargin )
%  EVAL - Evaluate Green function elements.
%
%  Usage for obj = stratified.tab.intra2 :
%    y = eval( obj, pos1, pos2, k0, varargin )
%  Input
%    pos1     :  observation points
%    pos2     :  source points
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    rows     :  positions of same size or not
%    smooth   :  w/o or with quasistatic contribution

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'smooth', 0 );
%  parse input
parse( p, varargin{ : } );

%  fill Green function tables ?
if isempty( obj.k0 ) || k0 ~= obj.k0,  obj = fill( obj, k0 );  end

%  convert to coordinates for stratified system
[ r, Z1, Z2, siz ] = fun( obj, pos1, pos2, varargin{ : } );
r( r < min( obj.rtab ) ) = min( obj.rtab );

%  perform interpolation
for name1 = [ "te", "tm" ]
for name2 = [ "", "z", "zz", "r", "rz", "s" ]
  y.( name1 + name2 ) = reshape(  ...
    obj.ytab1.( name1 + name2 )( r, Z1 ) +  ...
    obj.ytab2.( name1 + name2 )( r, Z2 ), siz );
end
end

%  add quasistatic contribution ?
if ~p.Results.smooth
  %  compute quasistatic Green function
  y2 = quasistatic( obj, pos1, pos2, k0, varargin{ : }, 'refl', 1 );
  %  add to tabulated values
  for name = convertCharsToStrings( fieldnames( y2 ) ) .'
    y.( name ) = y.( name ) + y2.( name ); 
  end
end


function [ r, Z1, Z2, siz ] = fun( obj,  pos1, pos2, varargin )
%  FUN - Convert to coordinates for stratified system.

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
    Z1 = bsxfun( @minus, pos1( :, 3 ), pos2( :, 3 ) .' ) + z2 - z1;
    Z2 = bsxfun( @plus,  pos1( :, 3 ), pos2( :, 3 ) .' ) - 2 * z1;
    %  size of output array
    siz = [ siz1( 1 : end - 1 ), siz2( 1 : end - 1 ) ];
  case 1
    %  radial distance and z-value
    r  = sqrt( ( pos1( :, 1 ) - pos2( :, 1 ) ) .^ 2 +  ...
               ( pos1( :, 2 ) - pos2( :, 2 ) ) .^ 2 );
    Z1 = pos1( :, 3 ) - pos2( :, 3 ) + z2 - z1;
    Z2 = pos1( :, 3 ) + pos2( :, 3 ) - 2 * z1;
    %  size of output array
    siz = siz1( 1 : end - 1 );
end
