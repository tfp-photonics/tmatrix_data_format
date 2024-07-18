path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;

wavelength = (400:50:800)';
epsilon = epsAu(wavelength); 
medium=1.33; % water

% spheroid semi-axes
stParams.a=20; stParams.c=40;

% simulation parameters
stParams.N=3; stParams.nNbTheta=120;

% additional options if required
stOptions = {};

% allocate 3D array for all results
qmax = 2*(stParams.N*(stParams.N + 1) + stParams.N );
tmatrix = zeros(length(wavelength), qmax, qmax);
% loop over wavelengths
for i=1:length(wavelength)
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;
    [stCoa, stT] = slvForT(stParams,stOptions);
    [T, u, up] = expand_tmat(stT, qmax); % tmatrix elements and indices
    ind = sub2ind(size(tmatrix),u*0+i, u, up); % linear index in 3D array
    tmatrix(ind) =  T(:,7) + 1i*T(:,8);
end
% keep note of which elements are known to be exact zeros
analytical_zeros = true(qmax,qmax);
analytical_zeros(sub2ind(size(analytical_zeros), u, up)) = false;

% data to export
epsilon = struct(...
    'embedding', medium^2, 'particle', epsilon, 'embedding_name', 'H2O, Water', ...
    'embedding_keywords', 'non-dispersive', 'embedding_reference', 'constant', ...
    'material_name', 'Au, Gold', 'material_reference', 'Au from Johnson and Christy 1972',...
    'material_keywords', 'dispersive, plasmonic');

geometry = struct('radiusxy', stParams.a, 'radiusz', stParams.c);

computation = struct('description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids', ...
    'accuracy', 1e-10, 'Lmax', stParams.N, 'Ntheta', stParams.nNbTheta, ...
    'analytical_zeros', analytical_zeros);

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'gold, spheroid, ebcm', 'script', [mfilename '.m']);

tmatrix_hdf5('smarties_simple.tmat.h5', tmatrix, wavelength, epsilon,...
              geometry, computation, comments)
