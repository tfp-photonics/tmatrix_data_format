%% stratified.layerstructure
%
% Base class for defining the layerstructure and for field computations.
%
%% Initialization

%  initialize layerstructure
%    mat    -  material vector in ascending order
%    z      -  vertical interface positions in ascending order
obj = stratified.layerstructure( mat, z );

%% Methods
%
% Once the layerstructure object has been initialized, the following
% methods are available:

%  number of interfaces
obj.n
%  locate positions in layer structure
ind = indlayer( obj, pos );
%  pairwise distance to layer scaled by bounding box radius
%    tau      -  boundary elememts
%    pos      -  positions
%    d        -  distance scaled by bounding box radius
%    iz       -  layer index for closest interface
[ d, iz ] = bdist2( obj, tau1, tau2 );
[ d, iz ] = bdist2( obj, tau1, pos2 );

%%
% We provide the following functions for the transfer matrix approach 

%  Fresnel coefficients between layers (i1,i2)
%    'efield'   -  TM coefficients for electric or magnetic (default) field
[ r, t ] = fresnel( obj, k0, kpar, i1, i2, PropertyPairs );
%  propagation matrix for given layer I1
p = propagate( obj, k0, kpar, i1 )
%  reflection and transmission coefficients for entire layer structure
%    'dir'    :  'up' for upgoing and 'down' for downgoing wave
[ r, t ] = rtcoeffs( obj, k0, kpar, PropertyPairs );
%  transfer matrix for given interface
m = transfer( obj, k0, kpar, i1, i2 );
%  total transfer matrix for layer structure
mtot = transfertot( obj, k0, kpar );

%%
% With these functions we can compute the electromagnetic fields within a
% layer structure through

%  secondary and primary fields 
%    k0         -  wavenumber of light in vacuum
%    kpar       -  parallel wavenumber
%    z1         -  field positions
%    z2         -  source position
%    'primary'  -  add primary field
f = fields( obj, k0, kpar, z1, z2, PropertyPairs );
%  propagate plane wave through layer structure
%    pol      -  light polarization
%    dir      -  light propagation direction
%    k0       -  wavenumber of light in vacuum
%    pos      -  positions where fields are evaluated
%   'primary' -  add primary fields ?
[ e, h ] = planewave( obj, pol, dir, k0, pos, PropertyPairs );
%  coefficients for secondary waves
%    i1     -  medium index for field point
%    i2     -  medium index for source point
wave = secondary( obj, k0, kpar, i1, i2 );

%% Grid functions
%
% The toolbox also provides a number of auxiliary functions for setting up
% tabulation grids for the reflected Green's functions. 

%  group position or position pairs into unique layer slices
pts = slice( obj, pos1 );
pts = slice( obj, pos1, pos2 );
%  tabulation grids for position pairs
%    'nr'   -  number of radial grid positons
%    'nz'   -  number of vertical grid positions
y = grid( obj, pts, PropertyPairs );

%%
% Consider for instance a substrate and a position vector extending over
% both layers.

%  glass-air substrate
layer = stratified.layerstructure( [ Material( 2.25, 1 ), Material( 1, 1 ) ], 0 ); 
%  positions
[ r, z ] = ndgrid( linspace( 0, 50, 51 ), linspace( -50, 50, 51 ) );
pos = [ r( : ), 0 * r( : ), z( : ) ];
%  group position in unique layer slices
pts = slice( layer, pos )

%%
%  pts = 
% 
%    1×2 struct array with fields:
% 
%     i1
%     pos
%     ind
%
% pts is a structure array with the layer index i1, the positions pos, and
% the index with respect to the calling function. Correspondingly, when
% slicing a position pair the positions are separated into 
% contributions within unique layers

%  group position pairs in unique layer slices
pts = slice( layer, pos, pos )

%%
%  pts = 
% 
%    1×4 struct array with fields:
% 
%     i1
%     i2
%     pos1
%     pos2
%     ind1
%     ind2
%
% Finally, we can pass the pts structure to the grid function in order to
% set up a grid that can be used for the tabulation of Green's functions.

%  tabulation grids for position pairs
r = grid( layer, pts )

%%
%  r = 
% 
%    1×4 struct array with fields:
% 
%     layer
%     i1
%     i2
%     r
%     z1
%     z2
%
%% Examples
%
% * <matlab:edit('demostratlayer01') demostratlayer01.m> |-| Primary and
% secondary fields in layerstructure.
% * <matlab:edit('demostratlayer02') demostratlayer02.m> |-|
% Electromagnetic fields for planewave excitation in layerstructure.
%
% Copyright 2023 Ulrich Hohenester

