function [ e, x ] = efield1( obj, far, varargin )
%  EFIELD - Electric field on image side using FFT.
%
%  Usage for obj = optics.lensimage2 :
%    [ e, x ] = efield1( obj, far, PropertyPairs )
%  Input
%    far    :  electric far-field along directions of obj.dir
%  PropertyName
%    n      :  size of output array
%    focus  :  focus position of imaging lens
%  Output
%     e     :  electric image fields in focal plane
%     x     :  image coordinates

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'n', size( obj.ind, 1 ) );
addParameter( p, 'focus', [ 0, 0, 0 ] );
%  parse input
parse( p, varargin{ : } );

%  refractive indices and wavenumbers
n1 = obj.mat1.n( obj.k0 );  k1 = obj.k0 * n1;
n2 = obj.mat2.n( obj.k0 );  k2 = obj.k0 * n2;
%  manipulate fields before reference sphere
far = far .* exp( 1i * k1 * obj.dir * p.Results.focus .' );
far = far * obj.rot;

dir = obj.dir * obj.rot;
[ up1, ut1, ur1 ] = optics.cart2unit( dir );
%  transformation for change of coordinate systems
far = optics.unitproject( far, ut1, ur1 ) +  ...
      optics.unitproject( far, up1, up1 );
%  apply manipulation function in backfocal plane
if ~isempty( obj.backfocal ), far = obj.backfocal( far );  end   

%  multiply with prefactor for Richards-Wolf integral
fac = 1i * k2 * ( n1 / n2 ) ^ 1.5 / ( 2 * pi ) ./ sqrt( dir( :, 3 ) );
far = far .* fac;
%  expand array to full size
[ i1, i2 ] = find( obj.ind );
n = size( obj.ind, 1 );
fun = @( k ) accumarray( { i1, i2 }, far( :, k ), [ n, n ] );
far = cat( 3, fun( 1 ), fun( 2 ), fun( 3 ) );

%  inverse Fourier transform
e = ifft2( far, p.Results.n, p.Results.n );
e = fftshift( fftshift( e, 1 ), 2 );
%  reshape and scale electric field
M = p.Results.n / n;
e = e * ( 2 * sin( obj.theta ) * M ) ^ 2;

%  image coordinates
x = 0.5 * pi * n * linspace( -1, 1, p.Results.n ) / ( k1 * sin( obj.theta ) );
%  apply shift phase
fac = exp( - 1i * k1 * sin( obj.theta ) * ( x + x .' ) );
e = reshape( fac( : ) .* reshape( e, [], 3 ), size( e ) );


