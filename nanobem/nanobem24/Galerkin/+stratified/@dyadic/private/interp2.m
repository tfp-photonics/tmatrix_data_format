function [ y1, y2 ] = interp2( tab, pos1, pos2, k0, varargin )
%  INTERP2 - Interpolate Green function table and compute derivatives.
%
%  Usage :
%    [ y1, y2 ] = interp2( tab, pos1, pos2, k0, PropertyPairs )
%  Input
%    tab      :  tabulated Green function
%    pos1     :  observation points, must be all located in same layer
%    pos2     :  source points, must be all located in same layer
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    eta      :  increment for finite differences
%  Output
%    y1       :  Green function values
%    y2       :  Green function derivatives

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'eta', 1e-2 );
%  parse input
parse( p, varargin{ : } );

%  increment vector 
eta = p.Results.eta;
incr = { [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] };
%  shift positions
Pos1 = cellfun( @( x ) bsxfun( @plus, pos1, eta * x ), incr, 'uniform', 0 );
Pos2 = cellfun( @( x ) bsxfun( @plus, pos2, eta * x ), incr, 'uniform', 0 );
%  assemble positions
Pos1 = vertcat( Pos1{ : }, pos1 );
Pos2 = vertcat( Pos2{ : }, pos2 );
%  interpolate Green function
z = interp( tab, Pos1, Pos2, k0, varargin{ : } );

for name = convertCharsToStrings( fieldnames( z ) ) .'
  %  reshape output
  zz = reshape( z.( name ), size( pos1, 1 ), 4, size( pos2, 1 ), 4 );
  %  values 
  z0 = zz( :, 4, :, 4 );
  %  first derivatives
  z1 = ( zz( :, 1 : 3, :, 4 ) - repmat( z0, 1, 3, 1, 1 ) ) / eta;
  z2 = ( zz( :, 4, :, 1 : 3 ) - repmat( z0, 1, 1, 1, 3 ) ) / eta;
  %  mixed derivative
  y2( 3 ).( name ) = ( zz( :, 1 : 3, :, 1 : 3 ) -  ...
               repmat( zz( :, 1 : 3, :,     4 ), 1, 1, 1, 3 ) -  ...
               repmat( zz( :,     4, :, 1 : 3 ), 1, 3, 1, 1 ) +  ...
               repmat( zz( :,     4, :,     4 ), 1, 3, 1, 3 ) ) / eta ^ 2;
  
  %  set output
  y1.( name ) = squeeze( z0 );
  %  first derivatives
  y2( 1 ).( name ) = permute( squeeze( z1 ), [ 1, 3, 2 ] );
  y2( 2 ).( name ) = squeeze( z2 );
end
