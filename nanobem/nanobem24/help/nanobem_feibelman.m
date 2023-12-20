%% Feibelman parameters
%
% Feibelman parameters allow for the incorporation of quantum effects into
% classical electrodynamic simulations. They were first introduced by
% Feibelman [1] for light reflection and transmission at planar interfaces,
% and were later reformulated by Yang et al. [2] in terms of modified
% boundary conditions. The main idea is to describe the fields on both
% sides of a boundary through solutions of the homogeneous Maxwell
% equations with local and isotropic permittivities, and to lump all
% modifications due to quantum effects at the boundary into two so-called
% Feibelman paramaters, which can be complex-valued and frequency
% dependent.
%
% Due to its formulation, the Feibelman approach is particularly well
% suited for the boundary element method approach. The nanobem toolbox
% provides a simple and versatile implementation of Feibelman parameters
% [3]. Simulations for single nanoparticles with fine discretizations using
% many boundary elements can be considerably slower than normal BEM
% simulations, say by a factor of three, but for coarse discretizations
% and/or for nanoparticles with separated boundaries the overhead of
% simulations including Feibelman parameters in comparison to normal BEM
% simulations is usually negligible.
%
% * [1] P. J. Feibelman, Prog. Surf. Sci. *12*, 287 (1982).
% * [2] Yi Yang et al., Nature *576*, 248 (2019).
% * [3] U. Hohenester and G. Unger, Phys. Rev. B *105*, 075428 (2022).
%
% Copyright 2022 Ulrich Hohenester
