%% stratified.dipole
%
% Excitation of oscillating dipole and enhancement of decay rate for
% oscillating dipole placed in stratified medium.
%
%% Initialization

%  initialization of oscillating dipole in stratified medium, additional
%  arguments will be passed to Sommerfeld integrator and tabulated Green
%  function object
%    layer      -  layer structure
%    pt         -  dipole positions, see stratified.Point
obj = stratified.dipole( layer, pt, PropertyPairs );

%% Methods

%  inhomogeneities for dipole excitation to be used for solution of BEM equations
%    tau  -  vector of boundary elements
%    k0   -  wavenumber of light in vacuum
qinc = eval( obj, tau, k0 );
qinc = obj( tau, k0 );
%  electromagnetic far-fields of oscillating dipole at Point positions PT
[ e, h ] = fields( obj, pt, k0 );
%  electromagnetic far-fields of oscillating dipole along directions DIR
[ e, h ] = farfields( obj, dir, k0 );
%  enhancement of total decay rate for solution SOL of BEM equations
%    tot    -  total decay rate
%    P0     -  decay rate of dipole in homogeneous medium
[ tot, P0 ] = decayrate( obj, sol );

%%
% The current toolbox functions are not optimized for accuracy and speed.
% The reflected dipole fields are obtained from the tabulated Green
% functions using finite differences, and the results are sometimes not
% overly accurate. In the decayrate function we compute the electric fields
% for all pairs of source and field positions, rather than just for the
% pairs of identical source and field positions, which can considerably
% slow down the simulations when many dipole positions are considered.
%
%% Examples
%
% * <matlab:edit('demostratdip01') demostratdip01.m> |-| Decay rate for
% dipole and gold sphere above substrate.
% * <matlab:edit('demostratdip02') demostratdip02.m> |-| LDOS enhancement
% for gold sphere above substrate.
%
% Copyright 2023 Ulrich Hohenester
