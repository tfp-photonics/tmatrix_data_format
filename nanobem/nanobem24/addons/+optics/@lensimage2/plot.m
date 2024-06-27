function plot( obj, far, varargin )
%  PLOT - Plot farfields on Gaussian reference sphere.
%
%  Usage for obj = lensimage :
%    plot( obj, far, PropertyPairs )
%  Input
%    far        :  electric far-fields along directions of obj.dir
%  PropertyName
%    colorbar   :  plot colorbar

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'colorbar', 0 );
%  parse input
parse( p, varargin{ : } );

%  expand array to full size
[ i1, i2 ] = find( obj.ind );
n = size( obj.ind, 1 );
fun = @( k ) accumarray( { i1, i2 }, far( :, k ), [ n, n ] );
far = cat( 3, fun( 1 ), fun( 2 ), fun( 3 ) );

%  functions
fun = { @real, @imag };
%  loop over Cartesian coordinates
for k = 1 : 3
for i = 1 : 2
  subtightplot( 2, 3, sub2ind( [ 3, 2 ], k, i ) );
  imagesc( fun{ i }( far( :, :, k ) ) .' );
  axis equal tight off
  title( [ 'k = ', num2str( k ), ', ', func2str( fun{ i } ) ] );
  if p.Results.colorbar,  colorbar;  end
end
end
