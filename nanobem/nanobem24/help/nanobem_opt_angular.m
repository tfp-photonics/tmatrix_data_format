%% optics.angular
%
% _optics.angular_ is a function that allows performing an angular spectrum
% decomposition, and is mainly implemented for testing purposes.
%
% * U. Hohenester, Nano and Quantum Optics (Springer, 2020), Eq. (3.8).
%
%% Initialization

%  initialize lens system for imaging
%     mat        -  material properties of embedding medium
%     k0         -  wavenumber of light in vacuum
%     dir        -  requested farfield directions
%     e          -  electric field E(x,y,z)
%     x,y        -  field positions
%     z=0        -  plane where field is computed
%    'nmax'      -  break up simulations into bunches of nmax=1000
far = angular( mat, k0, dir, e, x, y, z, 'nmax', nmax );

%%
% As output, optics.angular gives the far-field components propagating into
% directions dir.
%
%% Examples
%
% * <matlab:edit('demofocus01') demofocus01.m> |-| Focus fields for an incoming Laguerre-Gauss beam.
%
%%
%
% Copyright 2024 Ulrich Hohenester
