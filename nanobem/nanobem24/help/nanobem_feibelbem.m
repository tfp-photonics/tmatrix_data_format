%% feibelman.bemsolver
%
% feibelman.bemsolver objects solve the BEM equations including the
% mesoscopic boundary conditions based on Feibelman parameters. Its user
% interface and performance is similar to the galerkin.bemsolver object,
% with the only exceptions that the Feibelman parameters must be provided
% and the solutions have different tangential electromagnetic fields at the
% particle inside and outside.

%% Initialization

%  initialize BEM solver
%    TAU    -  boundary elements
%    PARAM  -  single feibelman.param object or array 
obj = feibelman.bemsolver( tau, param, PropertyPairs );

%%
% The additional parameters that can be passed in the form of property
% names and values to the BEM solver are discussed in more detail for the
% galerkin.bemsolver object.
%
%% Methods

%  solution of BEM equations with Feibelman parameters
%    QINC   - incoming fields, see e.g. galerkin.planewave
sol = solve( obj, qinc );
%  shorthand version for BEM solution
sol = obj \ qinc;
%  solution of BEM equations with and w/o Feibelman parameters
%    SOL1   -  solutiom with Feibelman parameters
%    SOL2   -  solution w/o  Feibelman parameters
[ sol1, sol2 ] = solve( obj, qinc );

%%
% |sol| and |sol1| are objects of type feibelman.solution, which hold the
% tangential electromagnetic fields at the particle inside and outside,
% |sol2| is a object of type galerkin.solution. The calling sequence with
% two solutions including and neglecting Feibelman parameters is often
% useful for the comparison of the different results, and is significantly
% faster than successive calls to feibelman.bemsolver and
% galerkin.bemsolver.
%
%% Examples
%
% * <matlab:edit('demofeibel01') demofeibel01.m> |-| Extinction cross
% section for nanosphere with constant Feibelman parameters.
% * <matlab:edit('demofeibel02') demofeibel02.m> |-| Field map for
% optically excited coupled spheres.
%
%%
%
% Copyright 2022 Ulrich Hohenester
