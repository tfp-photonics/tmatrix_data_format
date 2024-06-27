%% multipole.tmatrix
%
% multipole.tmatrix stores the T-matrices. For the multipoles we use the
% definition of Hohenester "Nano and Quantum Optics" (Springer, 2020), Eq.
% (E.4).
%
%% Initialization
%
% The T-matrix solver is initialized with

%  initialize solver for T-matrices
%    solver     -  T-matrix or Mie solver
%    k0         -  wavenumber of light in vacuum
tmat = multipole.tmatrix( solver, k0 );
%  add T-matrix components
tmat.aa = aa;     %  TM-TM elements
tmat.ab = ab;     %  TM-TE elements
tmat.ba = ba;     %  TE-TM elements
tmat.bb = bb;     %  TE-TE elements

%%
% Alternatively one can read the T-matrices stored in a h5-file, with the
% storage convention of  <https://github.com/tfp-photonics/treams treams>.

%  load T-matrix solver and T-matrices from input file
%    finp     -  h5 input file 
[ tsolver, tmat ] = h5load( finp );

%% Methods

%  convert between treams and nanobem conventions, used internally
tmat = convert( tmat, 'to_h5' );
tmat = convert( tmat, 'from_h5' );
%  expand T-matrices to full size [aa,ab;ba,bb]
t = full( tmat );
%  select specific T-matrix elements
%    fun    -  function @(l,m) for selection of specific (l,m) values
%    ind    -  index to selected elements of tmat.tab
tmat = select( tmat, fun );
tmat = select( tmat, ind );
%  average optical cross section
csca = scattering( tmat );
cext = extinction( tmat );

%%
% T-matrices can be used for computing the multipoles for a given incoming
% excitation.

%  inhomogeneity for incoming fields
%    fun        -  incoming fields [e,h]=fun(pos,k0) 
%    diameter   -  diameter for multipole expansion, Jackson Eq. (9.123)
q = qinc( tmat, fun );
q = qinc( tmat, fun, diameter );
%  solve T-matrix equation
%    sol        -  solution vector of type multipole.solution
sol = solve( tmat, q );

%% Save T-matrices
%
% T-matrices can be stored to a h5-file using the storage convention of
% <https://github.com/tfp-photonics/treams treams>. A more detailed
% description will be published in a forthcoming article by N. Asadova et
% al., _Suggestion for a T-matrix data format_, in preparation (2024).

%  save T-matrix or array of T-matrices to h5-file
%    fout     -  output file
%    info     -  additional info
h5save( tmat, fout, info );

%%
% The structure with the additional information can be produced with the
% multipole.h5info function.

%  obtain default info structure
info = multipole.h5info

%%
%  info = 
% 
%    struct with fields:
% 
%         description: ""
%               files: []
%            keywords: ""
%      matdescription: []
%             matname: []
%                name: ""
%     scatdescription: []
%       scatgroupname: []
%            scatname: []
%                 tau: []
%
% The entries can be modified by editing the fields of the structure or by
% passing their values in the form of property pairs to the call of
% multipole.hinfo. As an example, consider a T-matrix simulation for a
% dielectric nanosphere.

%  global name and description
info.name = "Dielectric nanosphere";
info.description = "Single sphere (160 nm) and single wavelength (1 mu)";
%  additional information for each material object
info.matname = [ "Embedding medium", "Dielectric material of sphere" ];
%  save discretized nanoparticle boundary in /computation/geometry
info.tau = tau;
%  save additional matlab scripts or data files in /computation, [] if none
info.files = [];

%%
% After creating the h5 data file _fout_ one might like to add additional
% information, such as the shape and size of basic three-dimensional
% geometries. This can be done by directly writing to the h5 file, as shown
% here for a dielectric sphere with a given diameter.

%  add shape and radius to scatterer (optional)
h5create( fout, '/scatterer_1/geometry/radius', 1 );
h5write( fout, '/scatterer_1/geometry/radius', 0.5 * diameter );
h5writeatt( fout, '/scatterer_1/geometry', 'shape', "sphere" );
h5writeatt( fout, '/scatterer_1/geometry/radius', 'unit', "nm" );

%% Examples
%
% * <matlab:edit('demomulti01') demomulti01.m> |-| T-matrix for dielectric
% nanosphere and single wavelength.
% * <matlab:edit('demomulti02') demomulti02.m> |-| T-matrices for TiO2
% nanodisk and multiple wavelengths.
%
% The following demo files show how to read in the previously stored
% T-matrices into a Python file to be used in combination with
% <https://github.com/tfp-photonics/treams treams>.
%
% * <matlab:edit('demomulti01.py') demomulti01.py> |-| Averaged scattering
% spectra for single nanoparticle.
% * <matlab:edit('demomulti02.py') demomulti02.py> |-| Same as
% demomulti02.py but for planewave excitation.
% * <matlab:edit('demomulti03.py') demomulti03.py> |-| Scattering spectra
% for coupled nanoparticles.
%
% Copyright 2024 Ulrich Hohenester

