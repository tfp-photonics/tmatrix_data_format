function h5save( obj, fout, varargin )
%  H5SAVE - Save T-matrices to file.
%
%  Usage for obj = multipole.tmatrix :
%    h5save( obj, fout, info )
%  Input
%    fout     :  output file
%    info     :  additional information, see tmatrix.h5info

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'info', multipole.h5info( varargin{ : } ) );
%  parse input
parse( p, varargin{ : } );

%  open output file
if isfile( fout ),  delete( sprintf( '%s', fout ) );  end
fid = H5F.create( fout );

%  write UUID to H5 file
id = matlab.lang.internal.uuid();
h5create( fout, '/uuid', 1, 'DataType', 'string' );
h5write( fout, '/uuid', id );
%  write name and description
info = p.Results.info;
for name = [ "name", "description" ]
  if info.( name ) ~= ""
    h5writeatt( fout, '/', name, info.( name ) );  
  end
end
%  additional information
h5writeatt( fout, '/', 'created_with',  ...
  "Matlab " + convertCharsToStrings( version( '-release' ) ) );
%  storage format version
[ v1, v2, v3 ] = H5.get_libversion();
ver = [ num2str( v1 ), '.', num2str( v2 ), '.', num2str( v3 ) ];
h5writeatt( fout, '/', 'storage_format_version', convertCharsToStrings( ver ) );

%  write wavenumber to H5 file
k0 = vertcat( obj.k0 );
h5create( fout, '/angular_vacuum_wavenumber', numel( k0 ) );
h5write( fout, '/angular_vacuum_wavenumber', k0 );
h5writeatt( fout, '/angular_vacuum_wavenumber', 'unit', "nm^{-1}" );
%  write T-matrix data
h5complex_write( fid, 'tmatrix', full( convert( obj, 'to_h5' ) ) );

%  permeabilities and permittivities
mat = obj( 1 ).solver.mat;
mu  = arrayfun( @( x ) x.mu ( k0 ), mat, 'uniform', 0 );
eps = arrayfun( @( x ) x.eps( k0 ), mat, 'uniform', 0 );
%  group name for materials
if info.matgroupname == ""
  info.matgroupname = arrayfun( @( i )  ...
    string( sprintf( 'Material%i', i ) ), 1 : numel( mat ), 'uniform', 1 );
end
%  write materials
plist = 'H5P_DEFAULT';
gid1 = H5G.create( fid, 'materials', plist, plist, plist );
for i = 1 : numel( mat )
  %  write permeability and permittivity to material group
  gid2 = H5G.create( gid1, info.matgroupname( i ), plist, plist, plist );
  h5complex_write( gid2, 'relative_permeability',  mu{ i } );
  h5complex_write( gid2, 'relative_permittivity', eps{ i } );
  close( gid2 );
  %  write attributes to material group
  name = "/materials/" + info.matgroupname( i );
  if ~isempty( info.matname )
    h5writeatt( fout, name, 'name', info.matname( i ) );
  end
  if ~isempty( info.matdescription )
    h5writeatt( fout, name, 'description', info.matdescription( i ) );
  end
end
%  close group identifier
close( gid1 );
%  create soft link to embedding medium
name = "/materials/" + info.matgroupname( obj( 1 ).solver.imat );
H5L.create_soft( name, fid, 'embedding', plist, plist );

%  table of spherical degrees and orders
tab = obj( 1 ).solver.tab;
%  write angular degrees 
n = numel( tab.l );
h5create( fout, '/modes/l', 2 * n, 'DataType', 'int64' );
h5write( fout, '/modes/l', int64( [ tab.l; tab.l ] ) );
%  write angular orders 
h5create( fout, '/modes/m', 2 * n, 'DataType', 'int64' );
h5write( fout, '/modes/m', int64( [ tab.m; tab.m ] ) );
%  write polarization
pol = repmat( [ "tm", "te" ], n, 1 );
h5create( fout, '/modes/polarization', 2 * n, 'DataType', 'string' );
h5write( fout, '/modes/polarization', pol( : ) );

%  write software
H5G.close( H5G.create( fid, 'computation', plist, plist, plist ) );
h5writeatt( fout, '/computation', 'software', "nanobem24" );
%  write contents of additional files
for i = 1 : numel( info.files )
  %  read file
  finp = fopen( info.files( i ) );
  str = convertCharsToStrings( fscanf( finp, '%c' ) );
  fclose( finp );
  %  write contents to H5 file
  name = sprintf( '/computation/file%i', i );
  h5create( fout, name, 1, 'DataType', 'string' );
  h5write( fout, name, str );
end
%  write geometry
if ~isempty( info.tau ),  h5geometry_write( fout, fid, info );  end
%  close file
H5F.close( fid );
