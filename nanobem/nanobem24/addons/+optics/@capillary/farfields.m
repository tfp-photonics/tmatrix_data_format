function [ e, h ] = farfields( obj, sol, dir, varargin )
%  FARFIELDS - Trace scattered light through capillary.  
%
%  Usage for obj = optics.capillary :
%    [ e, h ] = farfields( obj, sol, dir, PropertyPairs )
%  Input
%    sol        :  Mie solution for inner medium
%    dir        :  light propagation directions in outermost medium
%  PropertyName
%    shift      :  shift origin of particle
%  Output
%    e,h        :  electromagnetic farfields in outermost medium

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'shift', [] );
%  parse input
parse( p, varargin{ : } );

%  material parameters of capillary
mat = obj.mat;
%  wavenumber in vacuum and in media
k0 = sol.k0;
[ k1, kn ] = deal( mat( 1 ).k( k0 ), mat( end ).k( k0 ) );
%  wavenumber in z-direction and radial wavenumber in innermost medium
kz = kn * dir( :, 3 );
kr = sqrt( k1 ^ 2 - kz .^ 2 );
%  polar angle in innermost medium
t = atan2( kr, kz );

%  convert Mie coefficients from spherical to polar coordinates
[ a, b, M ] = sph2pol( sol.tab, sol.a, sol.b, t );
%  shift origin of particle
if ~isempty( p.Results.shift )
  [ a, b, M ] = shiftpol( a, b, M, t, k1 * p.Results.shift );
end
%  propagate coefficients from inner to outer medium
[ a, b ] = propagate( obj, a, b, k0, kz );

%  angles and directions of farfields
[ u, t ] = cart2sph( dir( :, 1 ), dir( :, 2 ), dir( :, 3 ) );
t = 0.5 * pi - t;
%  polar and azimuthal unit vectors
e1 = { cos( t ) .* cos( u ), cos( t ) .* sin( u ), - sin( t ) };
e2 = { - sin( u ), cos( u ), 0 * u };

%  prefactor
fac = bsxfun( @times, ( - 1i ) .^ ( M + 1 ), exp( 1i * M * u .' ) );
%  multiply coefficients with prefactor
a = reshape( bsxfun( @times, fac( : ), reshape( a, numel( fac ), [] ) ), size( a ) );
b = reshape( bsxfun( @times, fac( : ), reshape( b, numel( fac ), [] ) ), size( b ) );
%  sum over angular orders
siz = size( a );
a = reshape( sum( a, 1 ), [ siz( 2 : end ), 1 ] );
b = reshape( sum( b, 1 ), [ siz( 2 : end ), 1 ] );

%  electromagnetic farfields
fun = @( x, y ) bsxfun( @times, x, y );
e = cellfun( @( e1, e2 ) fun( b, e2 ) + fun( a, e1 ), e1, e2, 'uniform', 0 );
h = cellfun( @( e1, e2 ) fun( a, e2 ) - fun( b, e1 ), e1, e2, 'uniform', 0 );
%  put together farfields
[ k, Z ] = deal( sol.mat.k( k0 ), sol.mat.Z( k0 ) );
e = permute( cat( 3, e{ : } ), [ 1, 3, 2 ] ) / k * Z;
h = permute( cat( 3, h{ : } ), [ 1, 3, 2 ] ) / k;
