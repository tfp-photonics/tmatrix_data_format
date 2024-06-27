%% multipole.tsolver
%
% multipole.tsolver allows to compute the T-matrices for nanophotonic
% scatterers. For the multipoles we use the definition of Hohenester "Nano
% and Quantum Optics" (Springer, 2020), Eq. (E.4).
%
%% Initialization
%
% The T-matrix is initialized with

%  initialize solver for T-matrices
%    mat      -  material vector
%    imat     -  index for embedding medium      
%    lmax     -  maximal degree for multipole expansion
%   'rules'   -  quadrature rules
tsolver = multipole.tsolver( mat, imat, lmax, PropertyPairs );

%%
% Alternatively one can also specify a structure with selected elements for
% 'l' and 'm'.

%  initialize T-matrix solver with structure TAB
tsolver = multipole.tsolver( mat, imat, tab, PropertyPairs );
%  alternatively replace table for T-matrix solver
tsolver.tab = tab;

%% Methods
%
% The T-matrix solver has to be used in combination with a BEM solver such
% as <nanobem_ret_bem.html galerkin.bemsolver> or <nanobem_stratbem.html
% stratified.bemsolver>.  To compute the T-matrices for a given wavenumber,
% we first compute the solution for incoming multipole fields and then
% perform a multipole expansion of the solution vector.

%  solve BEM equation for multipole fields
%    bem    -  BEM solver
%    tau    -  boundary elements
%    k0     -  wavenumber of light in vacuum
sol = bem \ tsolver( tau, k0 );
%  perform multipole expansion of BEM solution
tmat = eval( tsolver, sol );

%%
% tmat is a T-matrix object of type multipole.solution that can be either
% stored for use in <https://github.com/tfp-photonics/treams treams> or
% used for various kinds of post-processing.
%
%% Examples
%
% * <matlab:edit('demomulti01') demomulti01.m> |-| T-matrix for dielectric
% nanosphere and single wavelength.
% * <matlab:edit('demomulti02') demomulti02.m> |-| T-matrices for TiO2
% nanodisk and multiple wavelengths.
%
% Copyright 2024 Ulrich Hohenester

