function h5geometry_write(fout, fid, gname,  info )
%  H5GEOMETRY_WRITE - Save geometry to H5 file.
%
%  Usage :
%    h5geometry_write( fout, fid, info )
%  Input
%    fout     :  output file
%    fid      :  file identifier for output file
%    gname    :  name of the scatterer group
%    info     :  additional information, see tmatrix.h5info

%  read auxiliary information
finp = fopen( 'h5readme.txt' );
str1 = convertCharsToStrings( fscanf( finp, '%c' ) );
fclose( finp );
H5G.create(fid, append(gname,'/geometry'), 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
H5G.create(fid, append(gname,'/geometry/mesh'), 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');

%  write auxiliary information to H5 file
space = H5S.create('H5S_SCALAR');  % Create a scalar dataspace
dtype = H5T.copy('H5T_C_S1');      % Copy the datatype for a string
H5T.set_size(dtype, 'H5T_VARIABLE'); % Set the datatype size to variable-length
% Create the dataset
H5D.create(fid, append(gname,'/geometry/mesh/readme'), dtype, space, 'H5P_DEFAULT');
%H5D.write(readme, dtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', {str});
h5write( fout, [append(gname,'/geometry/mesh/readme')], string(str1) );
%savegmsh( info.tau, 'mesh.msh' )
%  particle boundary and vertices
tau = info.tau;
%  unique vertices and faces
[ verts, ~, faces ] = unique( vertcat( tau.verts ), 'rows' );
faces = reshape( faces, [], 3 );
%  material index at boundary inside and outside
inout = vertcat( tau.inout );

%  write vertices
h5create( fout, append(gname,'/geometry/mesh/vertices'), size( verts ) );
h5write( fout, append(gname,'/geometry/mesh/vertices'), verts );
%  write faces
h5create( fout, append(gname,'/geometry/mesh/faces'), size( faces ) );
h5write( fout, append(gname,'/geometry/mesh/faces'), int64( faces ) );
%  write inside-outside informtaion
h5create( fout, append(gname,'/geometry/mesh/inout'), size( inout ) );
h5write( fout, append(gname,'/geometry/mesh/inout'), int64( inout ) );

% Accepted geometries
GEOMETRY_PARAMS.sphere = {'radius'};
GEOMETRY_PARAMS.ellipsoid = {'radiusx', 'radiusy', 'radiusz'};
GEOMETRY_PARAMS.spheroid = {'radiusxy', 'radiusz'};
GEOMETRY_PARAMS.cylinder = {'radius', 'height'};
GEOMETRY_PARAMS.cone = {'radius', 'height'};
GEOMETRY_PARAMS.torus = {'radius_major', 'radius_minor'};
GEOMETRY_PARAMS.helix = {'radius_helix', 'radius_wire', 'number_turns', 'pitch'};
GEOMETRY_PARAMS.cube = {'length'};
GEOMETRY_PARAMS.rectangular_cuboid = {'lengthx', 'lengthy', 'lengthz'};

if isfield(GEOMETRY_PARAMS, info.shape)
    shape_params = GEOMETRY_PARAMS.(info.shape);
    if length(shape_params) ~= length(fieldnames(info.params))
        error('Number of parameters provided does not match the expected number for shape %s.', info.shape);
    end

    for i = 1:numel(shape_params)
        fieldName = shape_params{i};
        fieldValue = info.params.(fieldName);
        if isnumeric(fieldValue)
            dataset = [append(gname,'/geometry/',fieldName)];
            h5writeatt( fout, append(gname,'/geometry'), 'shape', 'sphere' );
            h5create(fout, dataset, size(fieldValue));
            h5write(fout, dataset, fieldValue);
            h5writeatt( fout, dataset, 'unit', char(info.unit) );
        end
    end

else
    fprintf('Shape ''%s'' not found in GEOMETRY_PARAMS structure.\n', info.shape);
end

