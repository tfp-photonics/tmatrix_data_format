function h5save( obj, fout, info )
%  H5SAVE - Save T-matrices to file.
%
%  Usage for obj = multipole.tmatrix :
%    h5save( obj, fout, info )
%  Input
%    fout     :  output file
%    info     :  additional information, see tmatrix.h5info
%  open output file
if isfile( fout ),  delete( sprintf( '%s', fout ) );  end
fid = H5F.create( fout );
% %  write name and description
for name = [ "name", "description", "keywords" ]
  if info.( name ) ~= ""
    h5writeatt( fout, '/', name, char(info.( name )) );  
  end
end

%  write wavenumber to H5 file
k0 = squeeze(vertcat( obj.k0 ));

if size(k0, 1)
    h5create( fout, '/angular_vacuum_wavenumber', size( k0 ) );
    h5write( fout, '/angular_vacuum_wavenumber', k0 );
    h5writeatt( fout, '/angular_vacuum_wavenumber', 'unit', char("nm^{-1}") );
end
%  permeabilities and permittivities
mat = obj( 1 ).solver.mat;
mu  = arrayfun( @( x ) x.mu ( k0 ), mat, 'uniform', 0 );
eps = arrayfun( @( x ) x.eps( k0 ), mat, 'uniform', 0 );

%  write materials
plist = 'H5P_DEFAULT';
gide = H5G.create( fid, 'embedding', plist, plist, plist );
for i = 1 : numel( mat )
  %  write permeability and permittivity to material group
  if i == 1  
    h5complex_write( gide, 'relative_permeability',  mu{ i } );
    h5complex_write( gide, 'relative_permittivity', eps{ i } );
  else
    if numel( mat ) == 2
        grname = '/scatterer';
        gid1 = H5G.create( fid, grname, plist, plist, plist );
    else
        grname = compose('/scatterer_%d', i);
        gid1 = H5G.create( fid, grname, plist, plist, plist );
    end
    gid2 = H5G.create(fid,  append(grname,'/material'), plist, plist, plist );
    h5complex_write( gid2, 'relative_permeability',  mu{ i } );
    h5complex_write( gid2, 'relative_permittivity', eps{ i } );
    
%     %  write attributes to material group
    if ~isempty( info.matname )
    h5writeatt( fout, append(grname,'/material'), 'name', char(info.matname( i )) );
    end
    if ~isempty( info.matdescription )
    h5writeatt( fout, append(grname,'/material'), 'description', char(info.matdescription( i )) );
    end
    H5G.close( gid2 );
    
  end
end
H5G.close(gide);
% %  write geometry
if ~isempty( info.tau ),  h5geometry_write( fout, fid, grname, info );  end
%  close group identifier

close( gid1 );
%  table of spherical degrees and orders
tab = obj( 1 ).solver.tab;
%  write angular degrees 
n = numel( tab.l );
m = int64( repelem(tab.m, 2) );
l = int64( repelem(tab.l, 2) );
h5create( fout, '/modes/l', 2 * n, 'DataType', 'int64' );
h5write( fout, '/modes/l', l );
%  write angular orders 
h5create( fout, '/modes/m', 2 * n, 'DataType', 'int64' );
h5write( fout, '/modes/m', m );
%  write polarization
pol = repmat( [ "electric", "magnetic" ], 1, n );
h5create( fout, '/modes/polarization', 2 * n, 'DataType', 'string' );
h5write( fout, '/modes/polarization', pol( : ) );

%  write T-matrix data
%  change ordering to electric/magnetic
obj = convert( obj, 'to_h5' );
te = find( pol == "te" | pol == "magnetic" );
tm = find( pol == "tm" | pol == "electric" );
[ ~, ind ] = ismember( [ l( te ), m( te ) ], [ l( tm ), m( tm ) ], 'rows' ); 
te = te( ind );
if ndims(full(obj)) == 2
    data = zeros(size(full(obj), 1), size(full(obj), 2), 1);
else
    data = zeros(size(full(obj)));
end
for i = 1 : numel( obj )
  data( tm, tm, i ) = obj( i ).aa ;
  data( tm, te, i ) = obj( i ).ab;
  data( te, tm, i ) = obj( i ).ba;
  data( te, te, i ) = obj( i ).bb;
end
h5complex_write( fid, 'tmatrix', data);

%  write software
H5G.create( fid, '/computation', plist, plist, plist);
[ v1, v2, v3 ] = H5.get_libversion();
ver = [ num2str( v1 ), '.', num2str( v2 ), '.', num2str( v3 ) ];
software = append("nanobem24, ", "Matlab=" ,convertCharsToStrings( version( '-release' )), ", hdf5=",char(ver));
h5writeatt( fout, '/computation', 'software', char(software) );
h5writeatt( fout, '/computation', 'method', char("BEM, Boundary Element Method"));

%write main script
stack = dbstack('-completenames');
if length(stack) < 2
    error('This script must be run from another script.');
end
sourceScriptPath = stack(2).file;
[~, sourceScriptName, ext] = fileparts(sourceScriptPath);

gfiles = '/computation/files';
fileID = fopen(sourceScriptPath, 'r');
if fileID == -1
    error('Could not open source script file.');
end
source_script = fread(fileID, '*char')'; 
source_script = string(source_script);
fclose(fileID);

gid = H5G.create(fid, gfiles, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
H5G.close(gid);
dname = append(gfiles, '/', sourceScriptName, ext);
h5create(fout, dname, size(source_script), 'Datatype', 'string');
h5write(fout, dname, source_script); 

%  write contents of additional files
for i = 1 : numel( info.files )
  %  read file
  finp = fopen( info.files( i ) );
  str = convertCharsToStrings( fscanf( finp, '%c' ) );
  fclose( finp );
  %  write contents to H5 file
  name = sprintf( '/computation/files/%i', i );
  h5create( fout, name, 1, 'DataType', 'string' );
  h5write( fout, name, str );
end

gmt = '/computation/method_parameters/';
gid = H5G.create(fid, gmt, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
fields = fieldnames(info.method_parameters);
for j = 1:numel(fields)
    fname = fields{j};
    fdata = info.method_parameters.(fname);
    dataset = [append(gmt, fname)];
    h5create(fout, dataset, size(fdata), 'Datatype', class(fdata));
    h5write(fout, dataset, fdata); 
end
H5G.close(gid);
% %  close file;
H5F.close( fid );