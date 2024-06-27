%% optics.lensimage2
%
% Same as _optics.lensimage_ but using FFT. In the FFT implementation the
% image coordinates cannot be chosen independently but are returned by the
% field functions.  Note the different calling sequences in
% _optics.lensimage/efield_ and _optics.lensimage2/efield_.
%
%% Initialization

%  initialize lens system for imaging
%     mat1       -  material properties on object side
%     mat2       -  material properties on image  side
%     k0         -  wavenumber of light in vacuum
%     NA         -  numerical aperture
%    'rot'       -  rotation matrix for optical axis
%    'backfocal' -  field manipulation in backfocal plane
%    'n'         -  number of wavenumbers per dimension
lens = optics.lensimage2( mat1, mat2, k0, NA, PropertyPairs );

%%
% rot is an optional rotation matrix that allows rotating the optical axis
% from the +z-direction into other directions. backfocal is a user-defined
% function that can be used to manipulate the electromagnetic fields in the
% backfocal plane, for instance to mimic the effect of a quarter-wave
% plate.
%
%% Methods
%
% Suppose that we have the solution sol of a BEM simulation for a particle
% at hand, which is either of type _galerkin.solution_ for a homogeneous
% medium or stratified.solution for a stratified medium. To image the
% light scattered by the particle, we have to compute the farfields in the
% direction of the lens object and submit them to the Richards-Wolf
% approach.

%  imaging of farfields
%    sol      -  BEM solution for homogeneous or stratified medium
%    focus    -  shift focus position wrt origin of BEM simulation
%    n        -  size of output array, determines image resolution
%    x        -  image coordinates w/o magnification
far = farfields( sol, lens.dir );  %  compute farfields in directions lens.dir
plot( lens, far );                 %  plot farfields on Gaussian reference sphere
[ isca, x ] = efield( lens, far, 'focus', focus, 'n', n );  

%%
% Similarly, the fields of a planewave object ref, which is of type
% _optics.decompose_, can be imaged through

%  imaging of planewave decomposition
[ iref, x ] = efield( lens, ref, 'focus', focus, 'n', n ); 

%% Examples
%
% * <matlab:edit('demooptics02') demooptics02.m> |-| Imaging of point dipoles using optics.lensimage.
% * <matlab:edit('demooptics04') demooptics04.m> |-| Imaging of point dipoles using optics.lensimage2.
%
%%
%
% Copyright 2024 Ulrich Hohenester
