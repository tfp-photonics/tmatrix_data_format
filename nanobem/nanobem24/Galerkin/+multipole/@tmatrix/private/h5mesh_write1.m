function h5mesh_write1( fout, fid, info, imat )
%  H5MESH_WRITE1 - Save mesh to H5 file in scatterers.
%
%  Usage :
%    h5mesh_write1( fout, fid, info, imat )
%  Input
%    fout     :  output file
%    fid      :  file identifier for output file
%    info     :  additional information, see tmatrix.h5info
%    imat     :  material at outside of scatterers

%  particle boundary and materials at in- and outside
tau = info.tau;
inout = vertcat( tau.inout );

%  loop over materials of scatterers
for i1 = find( cellfun( @( x ) ~isempty( x ), imat, 'uniform', 1 ) )
  %  name of scatterer
  name = convertStringsToChars( "/" + info.matgroupname( i1 ) );
  %  vertices and faces
  ind = inout( :, 1 ) == i1;
  [ verts, ~, faces ] = unique( vertcat( tau( ind ).verts ), 'rows' );
  faces = reshape( faces, [], 3 );
  %  write vertices
  h5create( fout, [ name, '/geometry/mesh/vertices' ], size( verts ) );
  h5write(  fout, [ name, '/geometry/mesh/vertices' ], verts );
  h5writeatt( fout, [ name, '/geometry/mesh' ], 'unit', "nm" );
  %  write faces
  h5create( fout, [ name, '/geometry/mesh/faces' ], size( faces ) );
  h5write(  fout, [ name, '/geometry/mesh/faces' ], int64( faces ) );
  %  write position of scatterer
  h5create( fout, [ name, '/geometry/position' ], 3 );
  h5write(  fout, [ name, '/geometry/position' ], [ 0, 0, 0 ] );
  
  %  create soft links to outer material
  switch imat{ i1 }
    case 1
      name1 = "/embedding";
    otherwise
      name1 = "/" + info.matgroupname( i1 ) + "/material";
  end
  name2 = [ name, '/geometry/mesh/material_out' ];
  H5L.create_soft( name1, fid, name2, 'H5P_DEFAULT', 'H5P_DEFAULT' ); 

end
