%% optics.capillary
%
% The propagation of optical farfields through a capillary can be computed
% with objects of the _optics.capillary_ class. 
%
% * U. Hohenester et al., _Imaging the scattered light of a nanoparticle
% through a cylindrical capillary_, Nanophotonics 13, 457 (2024).
%
%% Initialization

%  initialization of capillary
%    mat        -  1×3 vector of material properties for capillary, from inside to outside
%    diameter   -  1×2 vector of diameters, from inside to outside
obj = optics.capillary( mat, diameter );

%%
% For a particle inside the capillary, whose optical properties are
% characterized by Mie coefficients, the farfields outside the capillary
% can be computed from

%  farfield propagation through capillary
%    sol    -  multipole.solution object with scattering coefficients a,b
%    dir    -  farfield propagation directions outside of capillary
%    shift  -  shift origin of multipole solution
[ e, h ] = farfields( obj, sol, dir );
[ e, h ] = farfields( obj, sol, dir, 'shift', shift );

%% Examples
%
% * <matlab:edit('democapillary01') democapillary01.m> |-| Side image of polystrene nanoparticle through capillary.
% * <matlab:edit('democapillary02') democapillary02.m> |-| Same as democapillary01.m but for different focus planes.
%
% Copyright 2024 Ulrich Hohenester
