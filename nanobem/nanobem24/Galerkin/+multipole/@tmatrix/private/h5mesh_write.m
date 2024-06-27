function h5mesh_write( fout, fid, info )
%  H5MESH_WRITE - Save mesh to H5 file.
%
%  Usage :
%    h5mesh_write( fout, fid, info )
%  Input
%    fout     :  output file
%    fid      :  file identifier for output file
%    info     :  additional information, see tmatrix.h5info

%  particle boundary and materials at in- and outside
tau = info.tau;
inout = vertcat( tau.inout );
%  determine for each inner material the outer materials
imat = accumarray( inout( :, 1 ), inout( :, 2 ), [], @( x ) { x } ) .';
imat = cellfun( @( x ) unique( x ), imat, 'uniform', 0 );

%  unique material combinations for boundaries ?
if all( cellfun( @( x ) numel( x ) <= 1, imat, 'uniform', 1 ) )
  %  save mesh information in scatterers
  h5mesh_write1( fout, fid, info, imat );
else
  %  save mesh information in /computation/mesh
  h5mesh_write2( fout, fid, info );
end
