function obj = init( obj, mat1, mat2, k0, NA, varargin )
%  INIT - Initialize imaging lens using FFT.
%
%  Usage for obj = optics.lensimage2 :
%    obj = init( obj, mat1, mat2, k0, NA, PropertyPair )
%  Input
%    mat1       :  material properties on object side
%    mat2       :  material properties on image  side
%    k0         :  wavenumber of light in vacuum
%    NA         :  numerical aperture
%  PropertyName
%    n          :  number of wavenumbers per direction
%    rot        :  rotation matrix for optical axis
%    backfocal  :  field manipulation in backfocal plane

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'n', 201 );
addParameter( p, 'rot', eye( 3 ) );
addParameter( p, 'backfocal', [] );
%  parse input
parse( p, varargin{ : } );

%  save input
[ obj.mat1, obj.mat2, obj.k0, obj.NA ] = deal( mat1, mat2, k0, NA );
%  rotation matrix and field manipulation in backfocal plane
obj.rot = p.Results.rot;
obj.backfocal = p.Results.backfocal;

%  opening angle of lens (check for consistency of NA)
n1 = mat1.n( k0 );
obj.theta = asin( NA / n1 );
assert( isreal( obj.theta ) );
%  grid for Gaussian sphere
x = sin( obj.theta ) * linspace( -1, 1, p.Results.n );
[ xx, yy ] = ndgrid( x );
dir = [ xx( : ), yy( : ), sqrt( 1 - xx( : ) .^ 2 - yy( : ) .^ 2 ) ];
%  keep only directions inside of cutoff angle
ind = true( size( xx ) );
ind( dir( :, 3 ) < cos( obj.theta ) ) = false;
%  save output
obj.dir = dir( ind, : ) * p.Results.rot .';
obj.ind = ind;
