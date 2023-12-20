%% multipole.miesolver
%
% Mie solver mainly for testing. The multipole.miesolver is initialized
% with

%  Initialize multipole.miesolver
%    mat1       -  material properties at sphere inside
%    mat2       -  material properties at sphere outside
%    diameter   -  diameter of sphere in nm    
%   'lmax'      -  maximum number of spherical degrees
mie = multipole.miesolver( mat1, mat2, diameter );
mie = multipole.miesolver( mat1, mat2, diameter, 'lmax', lmax );

%%
% We can then compute the T-matrix for the nanosphere with

%  T-matrix for nanosphere
%    k0       -  wavenumber of light in vacuum
%    tmat     -  multipole.tmatrix object
tmat = eval( mie, k0 );

%%
%
% Copyright 2023 Ulrich Hohenester

