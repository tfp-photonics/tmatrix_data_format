%% feibelman.solution
%
% The solutions of the BEM solver including Feibelman parameters are stored
% in the form of feibelman.solution objects. These objects can be used in
% the same way as galerkin.solution vectors.
%
%% Initialization
%
% The solution object is typically initialized in the call to the BEM
% solver including Feibelman parameters.

%  BEM is an object of type feibelman.bemsolver
%  QINC is a structure returned from an excitation class, 
%    e.g. galerkin.planewave or galerkin.dipole
sol = bem \ qinc;

%%
%  solution with properties:
%
%    inout: 2
%      tau: [1×284 BoundaryEdge]
%       k0: 0.0079
%        e: [426×1 double]
%        h: [426×1 double]
%
% In comparison to galerkin.solution objects, |sol| has the additional
% property |inout|, where 1 and 2 indicate that |e|, |h| are the fields at
% the particle inside or outside. To toggle between the inner and outer
% fields one can use

%  E,H are fields at particle inside
sol = toggle( sol, 1 );
%  E,H are fields at particle outside
sol = toggle( sol, 2 );

%%
% In general, the toggle function is only needed for the internal toolbox
% functions. Similarly to galerkin.solution objects, the Feibelman
% solutions can be used for excitation classes such as galerkin.planewave
% and galerkin.dipole, and for the computation of electromagnetic fields on
% and off the boundary.

%  electromagnetic fields at points PT1
%    'relcutoff'  -  cutoff for refined integration
%    'memax'      -  slice integration into bunches of size MEMAX
%    'waitbar'    -  show waitbar during evaluation
[ e, h ] = fields( sol, pt1, PropertyPairs );
%  evaluate tangential fields at centroid positions
[ e, h ] = interp( sol );
%  evaluate fields at boundary inside or outside
[ e, h ] = interp( sol, 'inout', inout );
%  surface charge distribution
sig = surfc( sol );

%% Examples
%
% * <matlab:edit('demofeibel01') demofeibel01.m> |-| Extinction cross
% section for nanosphere with constant Feibelman parameters.
% * <matlab:edit('demofeibel02') demofeibel02.m> |-| Field map for
% optically excited coupled spheres.
%
% Copyright 2022 Ulrich Hohenester
