function [ e, h ] = paraxial( obj, pos, k0, varargin )
%  PARAXIAL - Electromagnetic fields in paraxial approximation.
%    See Pan Song et al., J. Quant. Spec. & Rad. Trans. 241, 106713 (2020).
%
%  Usage for obj = laguerregauss :
%    [ e, h ] = paraxial( obj, pos, k0, PropertyPairs )
%  Input
%    pos    :  positions where fields are evaluated
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    shift  :  additional shift of positions
%  Output
%    e,h    :  electromagnetic fields in paraxial approximation

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'shift', [] );
%  parse input
parse( p, varargin{ : } );

%  additional shift of positions
if ~isempty( p.Results.shift )
  pos = tensor( pos, [ 1, 2 ] ) + tensor( p.Results.shift, [ 3, 2 ] );
  pos = reshape( double( pos, [ 1, 3, 2 ] ), [], 3 );
end

%  convert positions to cylinder coordinates
[ x, y ] = deal( pos( :, 1 ), pos( :, 2 ) );
[ u, r, z ] = cart2pol( x, y, pos( :, 3 ) );
r( r < 1e-10 ) = 1e-10;
%  wavenmumber and impedance
[ k, Z ] = deal( obj.mat.k( k0 ), obj.mat.Z( k0 ) );
%  azimuthal and radial index 
[ m, n ] = deal( obj.mLG, obj.nLG );

%  beam waist and Rayleigh range
w0 = waist( obj, k0 );
zr = 0.5 * k * w0 ^ 2;
w = w0 * sqrt( 1 + ( z / zr ) .^ 2 );
%  Laguerre polynomials
L1 = laguerre( n, m,     2 * r .^ 2 ./ w .^ 2 );
L2 = laguerre( n, m + 1, 2 * r .^ 2 ./ w .^ 2 );

%  amplitude, correct for missing normalization in obj.fun
a = 1i * amplitude( obj, k0 ) / ( 1i ^ m * 0.5 * ( w0 * k ) ^ 2 ); 
%  Eq. (3,4)
um = ( sqrt( 2 ) * r ./ w ) .^ m .*  ...
  exp( 1i * ( m * u - ( 2 * n + m ) * atan( z / zr ) ) );
u0 = 1 ./ ( 1 + 1i * z / zr ) .* exp( - ( r / w0 ) .^ 2 ./ ( 1 + 1i * z / zr ) );
ee = um .* u0 .* exp( 1i * k * z );

%  electric field, Eq. (37)
ez1 = ( m * ( 1i * x + y ) ./ ( k * r .^ 2 ) .* L1 -  ...
       1i * x ./ ( 1i * z - zr ) .* L1 - 4i * x ./ ( k * w .^ 2 ) .* L2 ) .* ee;
ez2 = ( m * ( 1i * y - x ) ./ ( k * r .^ 2 ) .* L1 -  ...
       1i * y ./ ( 1i * z - zr ) .* L1 - 4i * y ./ ( k * w .^ 2 ) .* L2 ) .* ee;     
%  magnetic field w/o impedance, Eq. (40)
hz1 = ( - m * ( 1i * x + y ) ./ ( k * r .^ 2 ) .* L1 +  ...
       1i * x ./ ( 1i * z - zr ) .* L1 + 4i * x ./ ( k * w .^ 2 ) .* L2 ) .* ee;
hz2 = ( + m * ( 1i * y - x ) ./ ( k * r .^ 2 ) .* L1 -  ...
       1i * y ./ ( 1i * z - zr ) .* L1 - 4i * y ./ ( k * w .^ 2 ) .* L2 ) .* ee;
     
%  polarization
[ polx, poly ] = polarization( obj, u );
%  assemble fields, Eqs. (35-40)
e = a * [   polx .* L1 .* ee, poly .* L1 .* ee, polx .* ez1 + poly .* ez2 ];
h = a * [ - poly .* L1 .* ee, polx .* L1 .* ee, poly .* hz1 + polx .* hz2 ] / Z;
%  deal with additional shift of positions
if ~isempty( p.Results.shift )
  e = permute( reshape( e, [], size( p.Results.shift, 1 ), 3 ), [ 1, 3, 2 ] );
  h = permute( reshape( h, [], size( p.Results.shift, 1 ), 3 ), [ 1, 3, 2 ] );
end
