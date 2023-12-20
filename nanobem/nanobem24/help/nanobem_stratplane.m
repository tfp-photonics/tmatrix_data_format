%% stratified.planewave
%
% Plane wave excitations and optical cross sections for optical plane wave
% excitations in stratified media are computed with the
% stratified.planewave class.
%
%% Initialization

%  initialize planewave object
%    layer      -  layerstructure
%    pol        -  polarization vectors
%    dir        -  light propagation directions
%   'rules'     -  quadrature rules 
obj = stratified.planewave( layer, pol, dir, PropertyPairs );

%% Methods
%
% After initialization, the stratified.planewave object can be used as
% follows.

%  inhomogeneity for BEM solver stratified.bemsolver
%    tau    -  vector of boundary elements
%    k0     -  wavenumber of light in vaccum
qinc = obj( tau, k0 );
qinc = eval( obj, tau, k0 );
%   electromagnetic fields at Point positions PT
[ e, h ] = fields( obj, pt, k0 );
%  absorption cross section for solution SOL of BEM equations
cabs = absorption( obj, sol );
%  extinction cross section for solution SOL of BEM equations
cext = extinction( obj, sol );
%  scattering cross section and differential cross section computed at
%  particle boundary SOL.TAU
[ csca, dsca ] = scattering( obj, sol );

%% Examples
%
% * <matlab:edit('demostratspec01') demostratspec01.m> |-| Spectrum for
% gold sphere above glass substrate.
% * <matlab:edit('demostratspec02') demostratspec02.m> |-| Spectrum for
% gold disk on top of glass substrate.
%
% Copyright 2023 Ulrich Hohenester
