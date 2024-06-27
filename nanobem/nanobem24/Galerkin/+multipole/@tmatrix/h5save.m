function h5save( obj, fout, info )
%  H5SAVE - Save T-matrices to file.
%
%  Usage for obj = multipole.tmatrix :
%    h5save( obj, fout, info )
%  Input
%    fout     :  output file
%    info     :  additional information, see tmatrix.h5info

%  open H5 output file
if isfile( fout ),  delete( sprintf( '%s', fout ) );  end
fid = H5F.create( fout );

%%
%  h5file
%   |
%   +-- name
%   +-- description
%   +-- keywords
%   +-- storage_format_version
%   +-- angular_vacuum_wavenumber

for name = [ "name", "description", "keywords" ]
  if info.( name ) ~= ""
    h5writeatt( fout, '/', name, info.( name ) );  
  end
end
h5writeatt( fout, '/', 'storage_format_version', "v1" ); 
%  write wavenumber to H5 file
k0 = vertcat( obj.k0 );
h5create( fout, '/angular_vacuum_wavenumber', numel( k0 ) );
h5write( fout, '/angular_vacuum_wavenumber', k0 );
h5writeatt( fout, '/angular_vacuum_wavenumber', 'unit', 'nm^{-1}' );

%%
%  h5file
%   |
%   +-- tmatrix
%   +-- modes
%        |
%        +-- l
%        +-- m
%        +-- polarization

%  combined indices for nanobem and H5 modes
ind1 = index( obj( 1 ).solver, [ 1, 2 ] );
ind2 = sortrows( index( multipole.base( max( ind1( :, 1 ) ) ), [ 1, 2 ] ) );
[ ~, i1 ] = ismember( ind2, ind1, 'rows' );
%  write T-matrix data
tmat = full( convert( obj, 'to_h5' ) );
h5complex_write( fid, 'tmatrix', tmat( i1, i1, : ) );
%  write angular degrees 
h5create( fout, '/modes/l', size( ind2, 1 ), 'DataType', 'int64' );
h5write( fout, '/modes/l', int64( ind2( :, 1 ) ) );
%  write angular orders 
h5create( fout, '/modes/m', size( ind2, 1 ), 'DataType', 'int64' );
h5write( fout, '/modes/m', int64( ind2( :, 2 ) ) );
%  write polarization
pol = [ "electric", "magnetic" ];
h5create( fout, '/modes/polarization', size( ind2, 1 ), 'DataType', 'string' );
h5write( fout, '/modes/polarization', pol( ind2( :, 3 ) ) );

%%
%  h5file
%    |
%    +-- embedding
%         |
%         +-- relative_permittivity
%         +-- relative_permeability
%         +-- name
%         +-- description
%    +-- NAME
%         |
%         +-- name
%         +-- description
%         +-- material
%              |
%              +-- relative_permittivity
%              +-- relative_permeability
%         +-- name
%         +-- description
%         +-- geometry

%  permeabilities and permittivities
mat = obj( 1 ).solver.mat;
mu  = arrayfun( @( x ) x.mu ( k0 ), mat, 'uniform', 0 );
eps = arrayfun( @( x ) x.eps( k0 ), mat, 'uniform', 0 );
%  material group names
if isempty( info.scatgroupname )
  info.matgroupname = [ "embedding",  ...
    arrayfun( @( i ) "scatterer_" + num2str( i ), 1 : numel( mat ) - 1, 'uniform', 1 ) ];
end

%  write embedding material properties
plist = 'H5P_DEFAULT';
gid = H5G.create( fid, '/embedding', plist, plist, plist );
h5complex_write( gid, 'relative_permeability',  mu{ 1 }, 'vector' );
h5complex_write( gid, 'relative_permittivity', eps{ 1 }, 'vector' );
%  write attributes to embedding group
if ~isempty( info.matname )
  h5writeatt( fout, '/embedding', 'name', info.matname( 1 ) );
end
if ~isempty( info.matdescription )
  h5writeatt( fout, '/embedding', 'description', info.matdescription( 1 ) );
end    
H5G.close( gid );

%  loop over remaining materials
for it = 2 : numel( mat )
  name = [ '/', convertStringsToChars( info.matgroupname( it ) ) ];
  
  H5G.create( fid, name, plist, plist, plist );
  %  write attributes to scatterer
  if ~isempty( info.scatname )
    h5writeatt( fout, name, 'name', info.scatname( it ) );
  end
  if ~isempty( info.scatdescription )
    h5writeatt( fout, name, 'description', info.scatdescription( it ) );
  end    
    
  %  write permeability and permittivity to material group
  gid = H5G.create(fid, [ name, '/material' ], plist, plist, plist );
  h5complex_write( gid, 'relative_permeability',  mu{ it }, 'vector' );
  h5complex_write( gid, 'relative_permittivity', eps{ it }, 'vector' );
  H5G.close( gid );
   
  %  write attributes to material group
  if ~isempty( info.matname )
    h5writeatt( fout, [ name, '/material' ], 'name', info.matname( it ) );
  end
  if ~isempty( info.matdescription )
    h5writeatt( fout, [ name, '/material' ], 'description', info.matdescription( it ) );
  end    
end

%%
%  h5file
%   |
%   +-- computation
%        |
%        +-- software
%        +-- method
%        +-- files

%  write software
H5G.create( fid, '/computation', plist, plist, plist );
[ v1, v2, v3 ] = H5.get_libversion;
ver = [ num2str( v1 ), '.', num2str( v2 ), '.', num2str( v3 ) ];
software = append( "nanobem24, ", "Matlab=" , ...
  convertCharsToStrings( version( '-release' ) ), ", hdf5=", ver );
h5writeatt( fout, '/computation', 'software', software );
h5writeatt( fout, '/computation', 'method', "BEM, Boundary Element Method" );

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

%  write mesh to scatterers or /computation
if ~isempty( info.tau ),  h5mesh_write( fout, fid, info );  end
%  close file
H5F.close( fid );
