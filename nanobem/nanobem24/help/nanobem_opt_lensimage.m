%% optics.lensimage
%
% _optics.lensimage_ allows to simulate the imaging of farfields and
% planewave decompositions, using the approach of Richards and Wolf. For
% detail see
%
% * U. Hohenester, Nano and Quantum Optics (Springer, 2020), Eq. (3.10).
%
%% Initialization

%  initialize lens system for imaging
%     mat1       -  material properties on object side
%     mat2       -  material properties on image  side
%     k0         -  wavenumber of light in vacuum
%     NA         -  numerical aperture
%    'rot'       -  rotation matrix for optical axis
%    'backfocal' -  field manipulation in backfocal plane
%    'nphi'      -  azimuthal discretization of Gaussian reference sphere
%    'ntheta'    -  polar discretization of Gaussian reference sphere
%    'mcut'      -  cutoff for angular degree
lens = optics.lensimage( mat1, mat2, k0, NA, PropertyPairs );

%%
% rot is an optional rotation matrix that allows rotating the optical axis
% from the +z-direction into other directions. backfocal is a user-defined
% function that can be used to manipulate the electromagnetic fields in the
% backfocal plane, for instance to mimic the effect of a quarter-wave plate
% as shown in the examples. nphi=51 and ntheta=50 are the discretiaztions
% of the Gaussian reference sphere used in the Richards-Wolf approach, and
% mcut allows for an optional truncation of the partial Fourier transform
% of the optical far-fields.
%
%% Methods
%
% Suppose that we have the solution sol of a BEM simulation for a particle
% at hand, which is either of type _galerkin.solution_ for a homogeneous
% medium or _stratified.solution_ for a stratified medium. To image the
% light scattered by the particle, we have to compute the farfields in the
% direction of the lens object and submit them to the Richards-Wolf
% approach.

%  imaging of farfields
%    sol      -  BEM solution for homogeneous or stratified medium
%    x        -  x-coordinates of image, w/o magnification
%    y        -  y-coordinates of image, w/o magnification
%    focus    -  shift focus position wrt origin of BEM simulation
far = farfields( sol, lens.dir );  %  compute farfields in directions lens.dir
plot( lens, far );                 %  plot farfields on Gaussian reference sphere
isca = efield( lens, far, x, y, 'focus', focus );  

%%
% Similarly, the fields of a planewave object ref, which is of type
% optics.decompose, can be imaged through

%  imaging of planewave decomposition
iref = efield( lens, ref, x, y, 'focus', focus ); 

%%
% isca and iref are arrays of size nx×ny×3.
%
%% Examples
%
% * <matlab:edit('demooptics02') demooptics02.m> |-| Imaging of point dipoles.
% * <matlab:edit('demooptics03') demooptics03.m> |-| Imaging of point dipoles, reversed optical axis.
% * <matlab:edit('demoiscat01') demoiscat01.m> |-| iSCAT images for a various nanoparticles.
% * <matlab:edit('demoiscat02') demoiscat02.m> |-| iSCAT images for gold nanosphere and varying focus planes.
% * <matlab:edit('demoiscat05') demoiscat05.m> |-| iSCAT images with interactive slider.
%
%%
%
% Copyright 2024 Ulrich Hohenester
