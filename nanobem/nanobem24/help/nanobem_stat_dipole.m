%% galerkinstat.dipole
%
% Excitation of oscillating dipole, and enhancement of total and radiative
% decay rate for oscillating dipole placed in photonic envrionment.
%
%% Initialization

%  initialization of oscillating dipole at Point position PT
%    pt           -  Point object for dipole positions
%   'relcutoff'   -  refined integration for close dipoles
%   'rules'       -  default quadrature rules
obj = galerkinstat.dipole( pt, PropertyPairs );

%% Methods

%  enhancement of total decay rate for given solution SOL of BEM equations
tot = decayrate( obj, sol );
%  inhomogeneities for dipole excitation to be used for solution of BEM equations
%    tau  -  vector of boundary elements
qinc = eval( obj, tau, k0 );
qinc = obj( tau, k0 );

%% Examples
%
% * <matlab:edit('demostatdip01') demostatdip01.m> |-| Lifetime of dipole
% above nanosphere.
%
% Copyright 2022 Ulrich Hohenester
