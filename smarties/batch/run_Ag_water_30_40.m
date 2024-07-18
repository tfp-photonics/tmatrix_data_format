path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;

%% Auto-generated file from _template.m

% requested precision (OA Cext)
accuracy = 1e-8;

% prolate Ag spheroid in water
% semi-axes a=b=20nm, c=40nm
% semi-axes a=b=30nm, c=40nm
wavelength = 350:5:850; wavelength = wavelength(:); 
Nl = length(wavelength);
epsilon=epsilon_Ag(wavelength);
medium=1.33;

% constant simulation parameters
stParams.a=30;
stParams.c=40;

% internal options
stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0; % NB will be estimated automatically
stOptions.bGetSymmetricT = false;
stOptions.bOutput = false; % verbosity

%% first, figure out the maximum N's needed

globalN = 1;
globalnNbTheta = 1;

% loop over wavelengths
for i=1:Nl
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;
    % Estimated convergence params
    [N, nNbTheta] = sphEstimateNandNT(stParams, stOptions, accuracy);
    stParams.N=N; stParams.nNbTheta=nNbTheta;
    % Increase params to test accuracy
    stParams2=stParams;
    stParams2.N=stParams2.N+5;
    stParams2.nNbTheta=stParams2.nNbTheta+5;

    [stCoa, stT] = slvForT(stParams,stOptions);
    [stCoa2, stT2] = slvForT(stParams2,stOptions);

    if(stOptions.bOutput)
        fprintf('Convergence testing... lambda = %.5g\n', wavelength(i));
        fprintf('<Cext> = %.10g,   relative error: %.2g\n', stCoa.Cext, abs(stCoa.Cext./stCoa2.Cext-1));
        fprintf('<Csca> = %.10g,   relative error: %.2g\n', stCoa.Csca, abs(stCoa.Csca./stCoa2.Csca-1));
        fprintf('<Cabs> = %.10g,   relative error: %.2g\n', stCoa.Cabs, abs(stCoa.Cabs./stCoa2.Cabs-1));
    end

    if(abs(stCoa.Cext./stCoa2.Cext-1) > 1.1*accuracy)
        % try once more

        stParams=stParams2;
        stParams.N=stParams2.N+5;
        stParams.nNbTheta=stParams2.nNbTheta+5;
        [stCoa, stT] = slvForT(stParams,stOptions);
        [stCoa2, stT2] = slvForT(stParams2,stOptions);

        if(abs(stCoa.Cext./stCoa2.Cext-1) > 1.1*accuracy)

            warning(fprintf('requested precision was not achieved for  = %.5g\n', wavelength(i)))
        end

    end

    globalN = max(globalN, stParams.N);
    globalnNbTheta = max(globalnNbTheta, stParams.nNbTheta);
end

%% now redo calculations for all wavelengths with these fixed params

stParams.N=globalN; stParams.nNbTheta=globalnNbTheta;
    
% allocate 3D array for all results
qmax = 2*(globalN*(globalN + 1) + globalN );
tmatrix = zeros(length(wavelength), qmax, qmax);

for i=1:Nl
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

%% data to export
epsilon = struct(...
    'embedding', medium^2, 'particle', epsilon, 'embedding_name', 'H2O, Water', ...
    'embedding_keywords', 'non-dispersive', 'embedding_reference', 'constant', ...
    'material_name', 'Ag', 'material_reference', 'Raschke et al 10.1103/PhysRevB.91.235137',...
    'material_keywords', 'dispersive');

geometry = struct('radiusxy', stParams.a, 'radiusz', stParams.c);

computation = struct('description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids', ...
    'accuracy', 1e-10, 'Lmax', stParams.N, 'Ntheta', stParams.nNbTheta, ...
    'analytical_zeros', analytical_zeros);

comments = struct('name', 'Ag prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'Ag, spheroid, ebcm', 'script', [mfilename '.m']);

tmatrix_hdf5('smarties_Ag_water_30_40.tmat.h5', tmatrix, wavelength, epsilon,...
              geometry, computation, comments)
