function obj = init( obj, mat1, mat2, k0, NA, varargin )
%  INIT - Initialize imaging lens.
%
%  Usage for obj = lensimage :
%    obj = init( obj, mat1, mat2, k0, NA, PropertyPair )
%  Input
%    mat1       :  material properties on object side
%    mat2       :  material properties on image  side
%    k0         :  wavenumber of light in vacuum
%    NA         :  numerical aperture
%  PropertyName
%    rot        :  rotation matrix for optical axis
%    backfocal  :  field manipulation in backfocal plane
%    nphi       :  number of azimuthal angles
%    ntheta     :  number of polar angles
%    mcut       :  cutoff for angular degree

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rot', eye( 3 ) );
addParameter( p, 'backfocal', [] );
addParameter( p, 'nphi',   51 );
addParameter( p, 'ntheta', 50 );
addParameter( p, 'mcut', inf );
%  parse input
parse( p, varargin{ : } );

%  save input
[ obj.mat1, obj.mat2, obj.k0, obj.NA ] = deal( mat1, mat2, k0, NA );
%  rotation matrix and field manipulation in backfocal plane
obj.rot = p.Results.rot;
obj.backfocal = p.Results.backfocal;

%  opening angle of lens, check for consistency of NA
n1 = mat1.n( k0 );
theta = asin( NA / n1 );
assert( isreal( theta ) );
%  angles for Gaussian sphere
n = 2 * fix( p.Results.nphi / 2 ) + 1;
obj.phi = ( 0 : ( n - 1 ) ) / n * 2 * pi;
[ obj.theta, obj.w ] = lgwt( p.Results.ntheta, 0, theta );
[ obj.theta, obj.w ] = deal( transpose( obj.theta ), transpose( obj.w ) );
%  cutoff for angular order
obj.mcut = min( p.Results.mcut, fix( p.Results.nphi / 2 ) );

%  directions for incoming far-fields
[ phi, theta ] = ndgrid( obj.phi, obj.theta );
obj.dir = cat( 2, cos( phi( : ) ) .* sin( theta( : ) ),  ...
                  sin( phi( : ) ) .* sin( theta( : ) ),  ...
                                     cos( theta( : ) ) ) * p.Results.rot .';
          