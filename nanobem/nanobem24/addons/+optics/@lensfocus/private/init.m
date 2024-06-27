function obj = init( obj, mat, k0, NA, varargin )
%  INIT - Initialize focusing lenses.
%
%  Usage for obj = lensfocus :
%    obj = init( obj, mat, NA, PropertyPair )
%  Input
%    mat      :  material properties on focal side
%    k0       :  wavenumber of light in vacuum
%    NA       :  numerical aperture
%  PropertyName
%    nphi     :  number of azimuthal angles
%    ntheta   :  number of polar angles
%    rot      :  rotation matrix for optical axis

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'nphi',   31 );
addParameter( p, 'ntheta', 30 );
addParameter( p, 'rot', eye( 3 ) );
%  parse input
parse( p, varargin{ : } );

%  save input
[ obj.mat, obj.k0, obj.NA ] = deal( mat, k0, NA );
%  rotation matrix
obj.rot = p.Results.rot;

%  opening angle of lens, check for consistency of NA
n = mat.n( k0 );
theta = asin( NA / n );
assert( isreal( theta ) );
obj.rad = sin( theta );
%  angles for fields after Gaussian sphere 
n = 2 * fix( p.Results.nphi / 2 ) + 1;
phi = ( 0 : ( n - 1 ) ) / n * 2 * pi;
[ theta, w ] = lgwt( p.Results.ntheta, 0, theta );
%  integration weights
w = ones( numel( phi ), 1 ) * diff( phi( 1 : 2 ) ) * w .';

%  grid of angles
[ phi, theta ] = ndgrid( phi, theta );
%  save angles and integration weights
[ obj.phi, obj.theta, obj.w ] = deal( phi( : ), theta( : ), w( : ) );
%  radii for system before Gaussian sphere
obj.rho = sin( obj.theta );
