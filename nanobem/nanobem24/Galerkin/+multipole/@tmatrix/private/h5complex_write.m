function h5complex_write( id, dataset, data )
%  H5COMPLEX_WRITE - Save complex matrix to dataset of H5 file.
%
%  Usage :
%    h5complex_write( id, dataset, data )
%  Input
%    id       :  group identifier
%    dataset  :  name of dataset
%    data     :  complex data for storage

%  convert value array to structure
wdata = struct( 'r', real( data ), 'i', imag( data ) );
%  create the complex datatype
doubleType = H5T.copy( 'H5T_NATIVE_DOUBLE' );
sz = H5T.get_size( doubleType );
%  compound datatype 
filetype = H5T.create ( 'H5T_COMPOUND', 2 * sz );
H5T.insert( filetype, 'r',  0, doubleType );
H5T.insert( filetype, 'i', sz, doubleType );

%  create the dataset and write the compound data to it
space = H5S.create_simple( ndims( data ), fliplr( size( data ) ), [] );
dset = H5D.create( id, dataset, filetype, space, 'H5P_DEFAULT' );
H5D.write( dset, filetype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', wdata );

%  close and release resources
H5D.close( dset );
H5S.close( space );
H5T.close( filetype );
