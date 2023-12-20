function h5geometry_write( fout, fid, info )
%  H5GEOMETRY_WRITE - Save geometry to H5 file.
%
%  Usage :
%    h5geometry_write( fout, fid, info )
%  Input
%    fout     :  output file
%    fid      :  file identifier for output file
%    info     :  additional information, see tmatrix.h5info

%  read auxiliary information
finp = fopen( 'h5readme.txt' );
str = convertCharsToStrings( fscanf( finp, '%c' ) );
fclose( finp );
%  write auxiliary information to H5 file
h5create( fout, '/geometry/readme', 1, 'DataType', 'string' );
h5write( fout, '/geometry/readme', str );

%  create soft links to materials
for i = 1 : numel( info.matgroupname )
  name1 = "/materials/" + info.matgroupname( i );
  name2 = sprintf( '/geometry/material%i', i );
  H5L.create_soft( name1, fid, name2, 'H5P_DEFAULT', 'H5P_DEFAULT' );
end

%  particle boundary and vertices
tau = info.tau;
%  unique vertices and faces
[ verts, ~, faces ] = unique( vertcat( tau.verts ), 'rows' );
faces = reshape( faces, [], 3 );
%  material index at boundary inside and outside
inout = vertcat( tau.inout );

%  write vertices
h5create( fout, '/geometry/vertices', size( verts ) );
h5write( fout, '/geometry/vertices', verts );
%  write faces
h5create( fout, '/geometry/faces', size( faces ) );
h5write( fout, '/geometry/faces', int64( faces ) );
%  write inside-outside informtaion
h5create( fout, '/geometry/inout', size( inout ) );
h5write( fout, '/geometry/inout', int64( inout ) );
