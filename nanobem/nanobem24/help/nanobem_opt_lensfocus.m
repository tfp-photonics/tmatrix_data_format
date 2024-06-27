%% optics.lensfocus
%
% _optics.lensfocus_ allows to simulate the focusing of incoming fields,
% using the approach of Richards and Wolf. For detail see
%
% * U. Hohenester, Nano and Quantum Optics (Springer, 2020), Eq. (3.13).
%
%% Initialization

%  initialize lens for focusing
%     mat        -  material properties on focus side
%     k0         -  wavenumber of light in vacuum
%     NA         -  numerical aperture
%    'rot'       -  rotation matrix for optical axis
%    'nphi'      -  azimuthal discretization of Gaussian reference sphere
%    'ntheta'    -  polar discretization of Gaussian reference sphere
lens = optics.lensfocus( mat, k0, NA, PropertyPairs );

%%
% Users can provide a rotation matrix rot that rotates the focused fields.
% The optional parameters nphi=31 and ntheta=30 control the discretization
% of the Gaussian reference sphere that is used to simulate focusing.
%
%% Methods
%
% After initialization, the user must provide the incoming fields impinging
% on the Gaussian reference sphere (with unit radius, the cutoff radius for
% a given NA is lens.rad) at the positions

%  requested positions for incoming fields, before crossing the Gaussian
%  reference sphere
[ x, y ] = deal( lens.x, lens.y );           %  Cartesian coordinates
[ phi, rho ] = deal( lens.phi, lens.rho );   %  polar coordinates

%%
% For instance, for an incoming Gauss-Laguerre beam with a topological
% charge of m=1 and x-polarization we use

%  electric field for Gauss-Laguerre beam
e = normpdf( lens.rho, 0, 0.5 ) .* exp( 1i * lens.phi );
e = e( : ) * [ 1, 0, 0 ];

%%
% Once the fields impinging on the Gaussian reference sphere are
% determined, the focused fields are obtained from

%  planewave decomposition of incoming fields
field = eval( lens, e );
field = eval( lens, e, 'focus', focus );    %  shift focus point

%%
% field is an optics.decompose object that can be evaluated for unbounded
% and stratified media, as discussed in the help pages for optics.decompose
% together with means to rotate and shift the focused fields.
%
%% Examples
%
% * <matlab:edit('demooptics01') demooptics01.m> |-| Focusing of incoming fields.
%
%%
%
% Copyright 2024 Ulrich Hohenester
