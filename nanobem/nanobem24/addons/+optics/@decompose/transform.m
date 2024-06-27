function obj = transform( obj, varargin )
%  TRANSFORM - Rotate or shift planewave decomposition.
%
%  Usage for obj = optics.decompose :
%    obj = transform( obj, 'rot', rot )
%    obj = transform( obj, 'shift', shift, 'mat', mat )
%  PropertyName
%    rot    :  rotation matrix
%    shift  :  shift vector
%    mat    :  material vector
%    
%  Output
%    obj    :  transformed planewave decomposition

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rot', [] );
addParameter( p, 'shift', [] );
addParameter( p, 'mat', [] );
%  parse input
parse( p, varargin{ : } );

%  rotate planewave decomposition
if ~isempty( p.Results.rot )
  obj.efield = obj.efield * p.Results.rot .';
  obj.dir = obj.dir * p.Results.rot .';
end
%  shift planewave decomposition
if ~isempty( p.Results.shift )
  k = p.Results.mat( obj.k0 );
  obj.efield = obj.efield .* exp( - 1i * k * obj.dir * p.Results.shift .' );
end
