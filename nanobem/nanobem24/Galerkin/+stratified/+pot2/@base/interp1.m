function [ y1, y2 ] = interp1( ~, green, pos1, pos2, k0, varargin )
%  INTERP1 - Interpolate tabulated Green function and derivatives.
%
%  Usage for obj = stratified.pot2.base :
%    [ y1, y2 ] = interp1( obj, green, pos1, pos2, k0, PropertyPairs )
%  Input
%    green  :  tabulated Green functions
%    pos1   :  observation points
%    pos2   :  source points
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    ind    :  tensor indices [i1,i2,q2,k]
%  Output
%    y1     :  reflected Green functions
%    y2     :  gradient of reflected Green function

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ind', 1 : 4 );
%  parse input
parse( p, varargin{ : } );

eta = 1e-4;
%  expand observation points
pos1 = cat( 3, bsxfun( @plus, pos1, eta * [ 1, 0, 0 ] ),  ...
               bsxfun( @plus, pos1, eta * [ 0, 1, 0 ] ),  ...
               bsxfun( @plus, pos1, eta * [ 0, 0, 1 ] ), pos1 );
pos1 = permute( pos1, [ 1, 3, 2 ] );
%  evaluate reflected Green functions
y = interp( green, pos1, pos2, k0, varargin{ : } );

if isempty( y )
  [ y1, y2 ] = deal( [] );
else
  %  dummy indices for tensor class
  ind = num2cell( p.Results.ind );
  [ i1, i2, q2, k ] = deal( ind{ : } );
  %  loop over reflected Green functions
  for name = convertCharsToStrings( fieldnames( y ) ) .'
    %  reflected Green function
    yy = y.( name );
    siz = size( yy );
    y1.( name ) = reshape( yy( :, 4, :, : ), siz( [ 1, 3, 4 ] ) );
    y1.( name ) = tensor( y1.( name ), [ i1, i2, q2 ] );
    %  derivative of reflected Green function
    y2.( name ) = cat( 2, yy( :, 1, :, : ) - yy( :, 4, :, : ),  ...
                          yy( :, 2, :, : ) - yy( :, 4, :, : ),  ...
                          yy( :, 3, :, : ) - yy( :, 4, :, : ) ) / eta;
    y2.( name ) = tensor( y2.( name ), [ i1, k, i2, q2 ] );                      
  end
end
