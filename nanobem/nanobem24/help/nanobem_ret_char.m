%% Characteristic modes
%
% Characteristic modes can be used for a simplified description of the
% optical response of nanoresonators.
%
% * Harrington & Mautz, IEEE Trans. Antennas and Prop. 5, 622 (1971).
% * Dey et al., ibid 68, 43 (2019).
%
%% Initialization
%
% To set up a solver for characteristic modes, one needs a BEM solver such
% as <nanobem_ret_bem.html galerkin.bemsolver> or <nanobem_stratbem.html
% stratified.bemsolver> to set up the solver

%  initialize solver for characteristic modes
%    bem    -  BEM solver galerkin.bemsolver or stratified.bemsolver
modes = charModes( bem );

%% Methods
%
% The characteristic modes for a given wavenumber can be computed and
% plotted with

modes = eval( modes, k0 );  %  evaluate characteristic modes
plot( modes );              %  plot characteristic modes
plot( modes, ind );         %  plot indexed characteristic modes

%%
% Once the characteristic modes are evaluated, we can solve the BEM
% equations using these modes or expand the solutions in terms of
% characteristic modes.

%  BEM solution using characteristic modes
%    q    -  inhomogeneities of incoming fields
%    ind  -  index for selected modes
sol = solve( modes, q );
sol = solve( modes, q, ind );
%  expansion coefficients for modes
a = expand( modes, sol );
a = expand( modes, sol, ind );
%  project BEM solution on characteristic modes
[ sol, a ] = project( modes, sol );
[ sol, a ] = project( modes, sol, ind );

%% Examples
%
% * <matlab:edit('demochar01') demochar01.m> |-| Characteristic modes for
% optically excited ellipsoid.
%
% Copyright 2023 Ulrich Hohenester

