%% stratified.tab.green
%
% Tabulated reflected Green's functions for matrix-friendly approach of
% Chew.
%
% * Chew et al., IEEE Antennas Wireless Propagt. Lett. 5 , 490 (2006).
% * Chew, Tong, and Hu, Integral equation methods for electromagnetic and
% elastic waves, Morgan and Claypool (2009).
%
%% Initialization
%
% To set up the tabulated reflected Green's function one first has to set
% up a layerstructure and define a computational grid using the toolbox
% function for gridding.

%  initialize tabulated reflected Green function object, the property pairs
%  are passed to the Sommerfeld integrator
%    grid       -  computational grid, see layerstructure/grid
y = stratified.tab.green( grid, PropertyPairs );

%%
% Consider a substrate and a nanosphere above the substrate. To set up the
% tabulated refletced Green function object, we proceed as follows

%  materials for glass susbtrate and air
mat1 = Material( 2.25, 1 );
mat2 = Material( 1, 1 );
%  layerstructure
layer = stratified.layerstructure( [ mat1, mat2 ], 0 );
%  nanosphere
p = trisphere( 144, 50 );
p.verts( :, 3 ) = p.verts( :, 3 ) - min( p.verts( :, 3 ) ) + 5;

%  slice positions and set up computational grid
r = slice( layer, p.pos, p.pos );
r = grid( layer, r, 'nr', 20, 'nz', 20 );
%  set up tabulated refletced Green function object
green = stratified.tab.green( r )

%%
%  green = 
% 
%    intra1 with properties:
% 
%         i1: 2
%       rtab: [1×20 double]
%       ztab: [1×20 double]
%      layer: [1×1 stratified.layerstructure]
%         k0: []
%     method: 'cubic'
%
% The function returns a tab.intra1 object that can be used for the
% tabulation of the reflected Green functions. Altogether, we provide three
% internal classes for tabulated Green functions.
%
% * *stratified.tab.intra1* both field and source points located in
% uppermost or lowest medium of layerstructure.
% * *stratified.tab.intra2* both field and source points located in
% the same medium of layerstructure.
% * *stratified.tab.inter* field and source points located in different
% media.
%
% Typically the computation of intra1 and intra2 objects is much faster
% that that of inter objects. In general green is a vector of different
% objects, corresponding to Green functions between field and source points
% in different layers.
%
%% Methods
%
% To fill the table of Green's functions one has to call

%  fill tabulated reflected Green functions
%    k0     -  wavenumber of light in vacuum
green = fill( green, k0 );

%%
% Depending on the number of Green's function objects this step can take
% from a few seconds to minutes. Once the Green's functions are filled, we
% can plot them

%  plot stratified.tab.intra1 object
plot( green );
%  plot stratified.tab.intra2 object
%    iz     -  plot table for (1) z1+z2 or (2) z1-z2
plot( green, iz );
%  plot stratified.tab.inter object
%    mode   -  fix (1) z2 or (2) z1
%    iz     -  fixed z2 or z1 value
plot( green, mode, 'iz', 1 );

%% 
% Alternatively, the tabulated Green function objects can be used for
% interpolation

%  interpolation of tabulated Green functions
%    pos1     -  field positions
%    pos2     -  source positions
%    k0       -  wavenumber of light in vacuum
y = interp( green, pos1, pos2, k0 );     
y = interp( green, pos1, pos2, k0, 'stat', 1 );    %  quasistatic contribution only
y = interp( green, pos1, pos2, k0, 'smooth', 1 );  %  smooth contribution only

%% Examples
%
% * <matlab:edit('demostratgreen01') demostratgreen01.m> |-| Tabulated
% Green function for sphere above substrate.
%
% Copyright 2023 Ulrich Hohenester

