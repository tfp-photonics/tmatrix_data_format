%% galerkin.dipole
%
% Excitation of oscillating dipole, and enhancement of total and radiative
% decay rate for oscillating dipole placed in photonic envrionment.
%
%% Initialization

%  initialization of oscillating dipole at Point position PT
%    pt           -  Point object for dipole positions
%   'relcutoff'   -  refined integration for dipoles located close to boundary
%   'rules'       -  default quadrature rules
%   'rules1'      -  integration rules for refinement
obj = galerkin.dipole( pt, PropertyPairs );

%% Methods

%  enhancement of total and radiative decay rate for given solution SOL of
%  BEM equations
[ tot, rad ] = decayrate( obj, sol );
%  inhomogeneities for dipole excitation to be used for solution of BEM equations
%    tau  -  vector of boundary elements
qinc = eval( obj, tau, k0 );
qinc = obj( tau, k0 );
%  electromagnetic far-fields of oscillating dipole along direction
%  specified by Point object PT
[ e, h ] = farfields( obj, pt, k0 );
%  electromagnetic far-fields of oscillating dipole at Point position PT
[ e, h ] = fields( obj, pt, k0 );

%% Examples
%
% * <matlab:edit('demoretdip01') demoretdip01.m> |-| Lifetime of dipole
% above nanosphere.
%
%
% Copyright 2022 Ulrich Hohenester

