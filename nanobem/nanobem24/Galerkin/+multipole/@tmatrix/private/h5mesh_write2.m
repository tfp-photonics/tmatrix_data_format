function h5mesh_write2( fout, fid, info )
%  H5MESH_WRITE - Save mesh to H5 file in /computation/mesh/.
%
%  Usage :
%    h5mesh_write2( fout, fid, info )
%  Input
%    fout     :  output file
%    fid      :  file identifier for output file
%    info     :  additional information, see tmatrix.h5info

%  read auxiliary information
finp = fopen( 'h5readme.txt' );
str = convertCharsToStrings( fscanf( finp, '%c' ) );
fclose( finp );
%  write auxiliary information to H5 file
h5create( fout, '/computation/mesh/readme', 1, 'DataType', 'string' );
h5write( fout, '/computation/mesh/readme', str );

%  particle boundary and material properties
tau = info.tau;
mat = tau( 1 ).mat;
%  create soft links to materials
for i = 1 : numel( mat )
  switch i
    case 1
      name1 = "/embedding";
    otherwise
      name1 = "/" + info.matgroupname( i ) + "/material";
  end
  name2 = sprintf( '/computation/mesh/material%i', i );
  H5L.create_soft( name1, fid, name2, 'H5P_DEFAULT', 'H5P_DEFAULT' );
end

%  unique vertices and faces
[ verts, ~, faces ] = unique( vertcat( tau.verts ), 'rows' );
faces = reshape( faces, [], 3 );
%  material index at boundary inside and outside
inout = vertcat( tau.inout );

%  write vertices
h5create( fout, '/computation/mesh/vertices', size( verts ) );
h5write( fout, '/computation/mesh/vertices', verts );
h5writeatt( fout, '/computation/mesh/vertices', 'unit', "nm" );
%  write faces
h5create( fout, '/computation/mesh/faces', size( faces ) );
h5write( fout, '/computation/mesh/faces', int64( faces ) );
%  write inside-outside informtaion
h5create( fout, '/computation/mesh/inout', size( inout ) );
h5write( fout, '/computation/mesh/inout', int64( inout ) );
