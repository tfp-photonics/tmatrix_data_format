%% stratified.bemsolver
%
% stratified.bemsolver is the central object for the BEM solution of the
% full Maxwell equations in presence of stratified media. It sets up the
% Calderon matrix and solves the BEM equations. Users not interested in the
% working details of the BEM solver can control its performance using just
% a few parameters, as described in <nanobem_ret_bem.html
% galerkin.bemsolver> and below.
%
%% Initialization

%  set up BEM solver for stratified media
%    tau          -  discretized particle boundary
%    layer        -  layer structure
obj = stratified.bemsolver( tau, layer, PropertyPairs );

%  the following property pairs can be passed to the BEM solver
%   'relcutoff'   -  relative cutoff for refined integration, see bdist2
%   'nduffy'      -  number of Legendre-Gauss points for Duffy integration
%   'rules'       -  default integration rules
%   'rules1'      -  integration rules for refinement
%   'order'       -  orders for series expansion, [] if none
%   'semi1'       -  scale real axis of semiellipse
%   'semi2'       -  scale imaginary axis of semiellipse 
%   'ratio'       -  z : r ratio which determines integration path
%   'op'          -  options for ODE integration 
%   'cutoff'      -  cutoff parameter for matrix-friendly approach of Chew
%   'nr'          -  number of radial tabulation points
%   'nz'          -  number of z-values for tabulation
%   'waitbar'     -  show waitbar during initialization

%%
% See <nanobem_ret_bem.html galerkin.bemsolver> for a description of the
% property pairs, as well as <nanobem_stratsomint.html
% stratified.isommerfeld> and <nanobem_stratpot1.html
% stratified.pot1.reflected>.
%
%% Methods
%
% Once the BEM solver is set up, we can compute the Calderon matrix or
% solve the BEM equations for some kind of external excitation.

%  fill matrices for tabulated Green function
%    k0     -  wavenumber of light in vacuum
obj = fill( obj, k0 );
%  compute Calderon matrix
cal = calderon( obj, k0 );
%  solve BEM equations for external excitation
%    q      -  structure with external excitation, e.g. stratified.planewave
[ sol, obj ] = solve( obj, q );
%  shortcut version 
sol = obj \ q;

%%
% The BEM solver returns an object of type <nanobem_stratsol.html
% stratified.solution> which can be used for various types of post
% processing.
%
%% Examples
%
% In the following we discuss the solution of a gold nanosphere above a
% glass-air substrate that is excited by a plane wave

%  materials for glass substrate, air, and gold
mat1 = Material( 2.25, 1 );
mat2 = Material( 1, 1 );
mat3 = Material( epstable( 'gold.dat' ), 1 );
%  layerstructure
layer = stratified.layerstructure( [ mat1, mat2 ], 0 );
%  nanosphere situated 5 nm above substrate
p = trisphere( 144, 50 );
p.verts( :, 3 ) = p.verts( :, 3 ) - min( p.verts( :, 3 ) ) + 5;
%  boundary elements with linear shape functions
tau = BoundaryEdge( [ mat1, mat2, mat3 ], p, [ 3, 2 ] );

%  set up BEM solver and planewave excitation
bem = stratified.bemsolver( tau, layer );
exc = stratified.planewave( layer, [ 1, 0, 0 ], [ 0, 0, -1 ] );
%  solve BEM equations 
k0 = 2 * pi / 600;
sol = bem \ exc( tau, k0 );

%  plot surface charges of nanopshere
plot( tau, surfc( sol ) );

%%
% <<../figures/bemstrat01.jpg>>
%
% * <matlab:edit('demostratspec01') demostratspec01.m> |-| Spectrum for
% gold sphere above glass substrate.
% * <matlab:edit('demostratspec02') demostratspec02.m> |-| Spectrum for
% gold disk on top of glass substrate.
%
% Copyright 2023 Ulrich Hohenester
