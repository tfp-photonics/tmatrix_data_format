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

%  convert between treams and nanobem conventions
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
sol = solve( tmat, qinc );

%% Save T-matrices
%
% T-matrices can be stored to a h5-file using the storage convention of
% <https://github.com/tfp-photonics/treams treams>.

%  save T-matrix or T-matrices to h5-file
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
%      matdescription: []
%        matgroupname: ""
%             matname: []
%                name: ""
%                 tau: []
%
% The entries can be modified by editing the structure of the fields or by
% passing their values in the form of property pairs to the call of
% multipole.hinfo. As an example, consider a T-matrix simulation for a
% dielectric nanosphere.

%  global name and description
info.name = "Dielectric nanosphere";
info.description = "Single sphere (160 nm) and single wavelength (1 mum)";
%  additional information for each material object
info.matgroupname = [ "embedding", "dielectric" ];
info.matname = [ "Embedding medium", "Dielectric material of sphere" ];
%  save discretized nanoparticle boundary in /geometry
info.tau = tau;
%  save additional matlab scripts or data files in /computation, [] if none
info.files = [];

%% Examples
%
% * <matlab:edit('demomulti01') demomulti01.m> |-| T-matrix for dielectric
% nanosphere and single wavelength.
% * <matlab:edit('demomulti02') demomulti02.m> |-| T-matrices for TiO2
% nanodisk and multiple wavelengths.
%
% Copyright 2023 Ulrich Hohenester

