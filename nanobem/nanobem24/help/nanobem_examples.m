%% Examples
%
% In the following we list the demo files provided by the nanobem toolbox,
% which are used elsewhere in the help pages.
%
%% Base classes
%
% * <matlab:edit('demomaterial01') demomaterial01.m> |-| Materials within nanobem toolobox.
% * <matlab:edit('demoparticle01') demoparticle01.m> |-| Plot elementary particle shapes.
% * <matlab:edit('demopoint01') demopoint01.m> |-| Place points in environment of sphere.
% * <matlab:edit('demopoint02') demopoint02.m> |-| Place points in environment of coupled spheres.
% * <matlab:edit('demoquad01') demoquad01.m> |-| Quadrature rules for triangle.
% * <matlab:edit('demoquadduffy01') demoquadduffy01.m> |-| Quadrature rules for identical and touching triangles.
% * <matlab:edit('demoquadduffy02') demoquadduffy02.m> |-| Apply quadrature rules for Duffy integration.
%
%% galerkin
%
% * <matlab:edit('demoboundaryedge01') demoboundaryedge01.m> |-| Plot Raviart-Thomas shape elements.
% * <matlab:edit('democalderon01') democalderon01.m> |-| Calderon matrix for nanosphere.
% * <matlab:edit('demoretdip01') demoretdip01.m> |-| Lifetime of dipole above nanosphere.
% * <matlab:edit('demoretfield01') demoretfield01.m> |-| Electric field of optically excited nanosphere.
% * <matlab:edit('demoretfield02') demoretfield02.m> |-| Electric field of optically excited coupled nanospheres.
% * <matlab:edit('demoretspec01') demoretspec01.m> |-| Optical spectrum for metallic nanosphere.
% * <matlab:edit('demoretspec02') demoretspec02.m> |-| Optical spectra for coupled nanospheres.
%
% _Simulations with resonance modes_
%
% * <matlab:edit('democim01') democim01.m> |-| Resonance modes for ellipsoid and planewave excitation.
% * <matlab:edit('democim02') democim02.m> |-| Resonance modes for ellipsoid and dipole excitation.
% * <matlab:edit('democim03') democim03.m> |-| Resonance modes for dielectric nanosphere and nonresonant background.
% * <matlab:edit('demochar01') demochar01.m> |-| Characteristic modes for optically excited ellipsoid.
%
% _Simulations with Feibelman parameters_
%
% * <matlab:edit('demofeibel01') demofeibel01.m> |-| Optical spectra nanospheres.
% * <matlab:edit('demofeibel02') demofeibel02.m> |-| Electric field of optically excited nanospheres.
% * <matlab:edit('demofeibel03') demofeibel03.m> |-| Resonance modes for ellipsoid and planewave excitation.
%
% _Simulations with stratified media_
%
% * <matlab:edit('demostratlayer01') demostratlayer01.m> |-| Primary and secondary fields in layerstructure.
% * <matlab:edit('demostratlayer02') demostratlayer02.m> |-| Electromagnetic fields for planewave excitation in layerstructure.
% * <matlab:edit('demostratsomint01') demostratsomint01.m> |-| Sommerfeld integration, plot integrand.
% * <matlab:edit('demostratsomint02') demostratsomint02.m> |-| Sommerfeld integration, plot integrand w/o quasistatic contribution.
% * <matlab:edit('demostratgreen01') demostratgreen01.m> |-| Tabulated Green function for sphere above substrate.
% * <matlab:edit('demostratspec01') demostratspec01.m> |-| Spectrum for gold sphere above glass substrate.
% * <matlab:edit('demostratspec02') demostratspec02.m> |-| Spectrum for gold disk on top of glass substrate.
% * <matlab:edit('demostratfield01') demostratfield01.m> |-| Electric field for optically excited nanosphere.
% * <matlab:edit('demostratfield02') demostratfield02.m> |-| Electric field for optically excited nanodisk.
% * <matlab:edit('demostratfar01') demostratfar01.m> |-| Electric farfields for optically excited nanosphere.
% * <matlab:edit('demostratfar02') demostratfar02.m> |-| Scattering spectra for detector with finite NA.
% * <matlab:edit('demostratdip01') demostratdip01.m> |-| Decay rate for dipole and gold sphere above substrate.
% * <matlab:edit('demostratdip02') demostratdip02.m> |-| LDOS enhancement for gold sphere above substrate.
%
% _Simulations with multipole expansion_
%
% * <matlab:edit('demomulti01') demomulti01.m> |-| T-matrix for dielectric nanosphere and single wavelength.
% * <matlab:edit('demomulti02') demomulti02.m> |-| T-matrices for TiO2 nanodisk and multiple wavelengths.
%
%% galerkinstat
%
% * <matlab:edit('demoboundaryvert01') demoboundaryvert01.m> |-| Plot linear shape elements.
% * <matlab:edit('demostatdip01') demostatdip01.m> |-| Lifetime of dipole above nanosphere.
% * <matlab:edit('demostatfield01') demostatfield01.m> |-| Electric field of optically excited nanosphere.
% * <matlab:edit('demostatfield02') demostatfield02.m> |-| Electric field of optically excited coupled nanospheres.
% * <matlab:edit('demostatspec01') demostatspec01.m> |-| Optical spectrum for metallic nanosphere.
% * <matlab:edit('demostatspec02') demostatspec02.m> |-| Optical spectra for coupled nanospheres.
%
% _Simulations with quasistatic eigenmodes_
%
% * <matlab:edit('demostateig01') demostateig01.m> |-| Quasistatic eigenmodes for nanosphere.
% * <matlab:edit('demostateig02') demostateig02.m> |-| Quasistatic eigenmodes for nanoellipsoid.
% * <matlab:edit('demostateig03') demostateig03.m> |-| Quasistatic eigenmodes for nanocube.
%
% Copyright 2022 Ulrich Hohenester