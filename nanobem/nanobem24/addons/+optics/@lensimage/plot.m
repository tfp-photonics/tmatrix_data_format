function varargout = plot( obj, far, varargin )
%  PLOT - Plot farfields on Gaussian reference sphere.
%
%  Usage for obj = lensimage :
%    h = plot( obj, far, PropertyPairs )
%  Input
%    far    :  electric far-fields along directions of obj.dir
%  Output
%     h     :  handle to plot object

%  sphere at infinity
[ u, t ] = ndgrid( [ obj.phi, 2 * pi ], obj.theta );
[ faces, verts ] =  ...
  surf2patch( sin( t ) .* cos( u ), sin( t ) .* sin( u ), cos( t ) );
pinfty = particle( verts * obj.rot, faces );
%  extend farfield
far = reshape( far, numel( obj.phi ), numel( obj.theta ), [] );
far( end + 1, :, : ) = far( 1, :, : );
far = reshape( far, numel( u ), [] );

%  plot farfields
h = plot( pinfty, far, varargin{ : } );
%  set output
if nargout == 1,  varargout{ 1 } = h;  end
