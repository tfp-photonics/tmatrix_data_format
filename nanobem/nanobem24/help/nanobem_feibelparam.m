%% feibelman.param
%
% feibelman.param objects define the Feibelman parameters for a set of
% boundary elements. Users must additionally provide a function that
% returns the perpendicular and parallel Feibelman parameters, using either
% a single value or multiple values for the different boundary elements.
%
%% Initialization

%  initialize Feibelman parameters
obj = feibelman.param( tau, fun );

%%
% |tau| is a vector of boundary elements defining the boundary for which
% the Feibelman parameters are defined. |fun| is a user-defined function
% that returns the perpendicular and parallel Feibelman parameters. For
% frequency independent parameters the function can be defined as

%  function for frequency-independent parameters
fun = @( ~ ) deal( dperp, dpar );

%%
% Otherwise the user must provide a function 
%  
%  function [ dperp, dpar ] = fun( ene )
%    %  FUN - User-defined function for Feibelman parameters.
%    %
%    %  Input
%    %    photon energy in eV
%    %  Output
%    %    dperp  :  perpendicular Feibelman parameter in nanometers
%    %    dpar   :  parallel Feibelman parameter in nanometers
%    ...
% end
%
%%
% The function can either return single values for |dperp| and |dpar|,
% which are used for all boundary elements |tau|, or vectors of the same
% length as |tau| containing different Feibelman parameters for the
% different boundary elements.
%
%% Methods
%
% The Feibelman parameters can be evaluated either for the set of boundary
% elements |obj.tau| or for a user-defined vector of boundary elements. In
% the latter case only boundary elements contained in |obj.tau| have
% Feibelman parameters with a non-zero value, all other Feibelman
% parameters are set to zero, such that the usual boundary conditions of
% tangential electromagnetic fields apply.

%  evaluate Feibelman parameters for a single OBJ
%    k0 is the wavenumber of light in vacuum
[ dperp, dpar ] = eval( obj, k0 );
%  evaluate Feibelman parameters for an array OBJ and for user-defined
%  boundary elements TAU
[ dperp, dpar ] = eval( obj, tau, k0 );

%% Example
%
% * <matlab:edit('demofeibel01') demofeibel01.m> |-| Extinction cross
% section for nanosphere with constant Feibelman parameters.
%
% Copyright 2022 Ulrich Hohenester

