function y = eval( obj, pos1, pos2, k0, varargin )
%  EVAL - Evaluate Green function elements.
%
%  Usage for obj = stratified.tab.inter :
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
[ r, z1, z2, siz ] = fun( obj, pos1, pos2, varargin{ : } );
r( r < min( obj.rtab ) ) = min( obj.rtab );

%  perform interpolation
for name = convertCharsToStrings( fieldnames( obj.ytab ) ) .'
  y.( name ) = reshape( obj.ytab.( name )( r, z1, z2 ), siz );
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


function [ r, z1, z2, siz ] = fun( ~, pos1, pos2, varargin )
%  FUN - Convert to coordinates for stratified system

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rows', 0 );
%  parse input
parse( p, varargin{ : } );

%  reshape position arrays
[ siz1, pos1 ] = deal( size( pos1 ), reshape( pos1, [], 3 ) );
[ siz2, pos2 ] = deal( size( pos2 ), reshape( pos2, [], 3 ) );

switch p.Results.rows
  case 0
    %  radial distance and z-values
    r = pdist2( pos1( :, 1 : 2 ), pos2( :, 1 : 2 ) );
    [ z1, z2 ] = ndgrid( pos1( :, 3 ), pos2( :, 3 ) );
    %  size of output array
    siz = [ siz1( 1 : end - 1 ), siz2( 1 : end - 1 ) ];
  case 1
    %  radial distance and z-values
    r = sqrt( ( pos1( :, 1 ) - pos2( :, 1 ) ) .^ 2 +  ...
              ( pos1( :, 2 ) - pos2( :, 2 ) ) .^ 2 );
    z1 = pos1( :, 3 );
    z2 = pos2( :, 3 );
    %  size of output array
    siz = siz1( 1 : end - 1 );
end
