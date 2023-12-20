%% feibelman.cimsolver
%
% The nanobem toolbox provides a CIM solver for resonance modes including
% Feibelman parameters.
%
%% Initialization
%
% For a given set of boundary elements |tau| and an array of Feibelman
% parameters |param|, we set up the CIM solver through

%  initialize CIM solver
cim = feibelman.cimsolver( tau, param, 'nr', 150 );
%  initialize CIM contour
obj.contour = cimbase.ellipse( [ zmin, zmax ], irad, n );

%% Methods
%
% The computation of the resonance modes can be done in the same way as
% described for the |galerkin.cimsolver| object. The only difference is
% that the CIM solver now only returns the fields at the particle outside.
% To additionally obtain the fields at the particle inside one has to
% proceed as follows.

%  solve BEM equations with Feibelman parameters using resonance modes
sol = cim \ qinc;
%  compute tangential electromagnetic fields at particle inside, 
%  only if requested, SOL is converted to a feibelman.solution object
sol = match( cim, sol );

%% Examples
%
% * <matlab:edit('demofeibel03') demofeibel03.m> |-| Resonance modes for
% nanoellipsoid using Feibelman parameters.
%
% Copyright 2022 Ulrich Hohenester
