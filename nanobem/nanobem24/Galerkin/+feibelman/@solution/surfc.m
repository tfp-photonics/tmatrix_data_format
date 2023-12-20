function sig = surfc( obj )
%  SURFC - Surface charge at centroids, Hohenester Eq. (7.24).
%
%  Usage for obj = feibelman.solution :
%    sig = surfc( obj )
%  Output
%    sig    :  surface charge

%  electric fields at centroids and particle inside and outside
e1 = interp( obj, 'inout', 1 );
e2 = interp( obj, 'inout', 2 );
%  dot product with outer surface normal
nvec = vertcat( obj.tau.nvec );
sig = bsxfun( @times, reshape( e2 - e1, numel( nvec ), [] ), nvec( : ) );
sig = squeeze( sum( reshape( sig, size( e1 ) ), 2 ) );
