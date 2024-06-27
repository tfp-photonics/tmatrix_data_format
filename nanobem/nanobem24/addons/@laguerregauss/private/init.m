function obj = init( obj, mat, varargin )
%  INIT - Initialize Laguerre-Gauss beam.
%
%  Usage for obj = laguerregauss :
%    obj = init( obj, mat, PropertyPairs )
%  Input
%    mat    :  material properties of embedding medium
%  PropertyName
%    foc    :  effective focal length
%    w0     :  input beam waist
%    pol    :  polarization before focusing
%    nLG    :  radial index
%    mLG    :  azimuthal index
%    pow    :  laser power (W)
%    M2     :  M-squared value     

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'foc', [] );
addParameter( p, 'w0',  [] );
addParameter( p, 'pol', [ 1, 0, 0 ] );
addParameter( p, 'nLG', 0 );
addParameter( p, 'mLG', 0 );
addParameter( p, 'pow', 1 );
addParameter( p, 'M2',  1 );
%  parse input
parse( p, varargin{ : } );

obj.mat = mat;
for name = { 'foc', 'w0', 'pol', 'nLG', 'mLG', 'pow', 'M2' }
  obj.( name{ 1 } ) = p.Results.( name{ 1 } );
end
