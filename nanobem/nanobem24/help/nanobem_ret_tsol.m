%% multipole.solution
%
% multipole.solution stores the multipole solution and allows for various
% kinds of post-processing.
%
%% Initialization
%
% The T-matrix solution is initialized with

%  Initialize T-matrix solution
%    lmax     -  maximal degree for multipole expansion
%    mat      -  material properties of embedding medium
%    k0       -  wavenumber of light in vacuum
%    a        -  TM scattering coefficients
%    b        -  TE scattering coefficients
%    ai       -  TM incoming coefficients
%    bi       -  TE incoming coefficients    
sol = multipole.solution( lmax, mat, k0, a, b, ai, bi );

%%
% Alternatively, one can initialize the object with fewer variables and add
% them later to the object.

%  Initialize T-matrix solution
sol = multipole.solution( lmax, mat, k0 );
sol.a = a;
sol.b = b;

%% Methods

%  electromagnetic multipole fields in embedding medium
%    pos    -  positions where fields are evaluated
[ e, h ] = fields( sol, pos );
%  electromagnetic far-fields
%    dir    -  propagation directions for far-fields
%    shift  -  shift origin of multipole expansion
[ e, h ] = farfields( sol, dir );
[ e, h ] = farfields( sol, dir, 'shift', shift );
%  scattering and extinction power, cext is divided by incoming power
csca = scattering( sol );
cext = extinction( sol );    
%  scattering power for user-defined detector
%    pinfty   -  face-vertex structure for detector
csca = scattering( sol, 'pinfty', pinfty );
%  optical force in pN and torque in pN×nm
[ f, n ] = optforce( sol );
%  re-use auxiliary data for multiple calls
[ f, n, data ] = optforce( sol );
[ f, n ] = optforce( sol. data );

%% Examples
%
% * <matlab:edit('demomulti02') demomulti02.m> |-| T-matrices for TiO2
% nanodisk and multiple wavelengths.
%
% Copyright 2024 Ulrich Hohenester

