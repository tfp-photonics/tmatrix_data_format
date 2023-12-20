%  DEMOFEIBEL03 - Compute resonance modes using Feibelman parameters, and
%  perform simulations for optically excited ellipsoid.

%  nanoellipsoid
diameter = 20;
p = trisphere( 144, diameter );
p.verts( :, 1 ) = p.verts( :, 1 ) * 2;

mat1 = Material( 1, 1 );
mat2 = Material( epsdrude( 'Au' ), 1 );
%  material properties
mat = [ mat1, mat2 ];

%  boundary elements with linear shape functions
tau = BoundaryEdge( mat, p, [ 2, 1 ] );
%  constant Feibelman parameters
[ d1, d2 ] = deal( - 0.4, 0 );
param = feibelman.param( tau, @( ~ ) deal( d1, d2 ) );

%  countour integral method solver
cim = feibelman.cimsolver( tau, param, 'nr', 150 );
cim.contour = cimbase.ellipse( [ 1.5, 2.6 ], 0.2, 60 );
%  compute contour integral and eigenvalues
rules = quadboundary.rules( 'quad3', triquad( 1 ) );
data = eval1( cim, 'rules', rules, 'waitbar', 1 );

tol = tolselect( cim, data, 'tol', 1e-2 );
if isempty( tol ),  return;  end
cim = eval2( cim, data, 'tol', tol );

% plot( horzcat( cim.contour.z ) );  hold on
% plot( cim.ene, 'm+' );

%%
%  planewave excitation
exc = galerkin.planewave( [ 1, 0, 0 ], [ 0, 0, 1 ] );

%  light wavelength in vacuum
lambda = linspace( 400, 700, 20 );
k0 = 2 * pi ./ lambda;
%  allocate optical cross sections
[ cext1, cext2 ] = deal( zeros( numel( k0 ), size( exc.pol, 1 ) ) );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  loop over wavenumbers
for i = 1 : numel( k0 )
  %  solution of BEM equations, fields at particle outside only
  sol1 = cim \ exc( tau, k0( i ) );
  cext1( i, : ) = extinction( exc, sol1 );
  
  sol2 = data.bem \ exc( tau, k0( i ) );
  cext2( i, : ) = extinction( exc, sol2 );
  
  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );

%%  final plot
figure
plot( lambda, cext1, 'o-' );  hold on
plot( lambda, cext2, '+-' );

legend( 'CIM', 'BEM' );

xlabel( 'Wavelength (nm)' );
ylabel( 'Cross section (nm^2)' );
