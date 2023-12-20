function [ csca, dsca ] = scattering( obj, sol )
%  SCATTERING - Scattering cross section.
%
%  Usage for obj = stratified.planewave :
%    [ csca, dsca ] = scattering( obj, sol )
%  Input
%    sol    :  solution of BEM equations
%  Output
%    csca   :  scattering cross section
%    dsca   :  differential cross section

%  allocate output for differential cross section
m = size( obj.pol, 1 );
dsca = zeros( numel( sol.tau ), m );

%  dummy indices for internal vector class
[ i, q, k, ipol ] = deal( 1, 2, 3, 4 );

%  loop over boundary elements
for pt = quadboundary( sol.tau, obj.rules )
  
  if pt.inout( 2 ) <= obj.layer.n + 1
    %  quadrature weights
    [ ~, w ] = eval( pt );
    w = tensor( w, [ i, q ] );
    %  normal vector
    nvec = tensor( vertcat( pt.tau.nvec ), [ i, k ] );
    %  interpolate electromagnetic fields
    [ e, h ] = interp( sol, pt );
    e = tensor( e, [ i, q, k, ipol ] );
    h = tensor( h, [ i, q, k, ipol ] );
    %  incoming electromagnetic fields
    [ einc, hinc ] = fields2( obj, pt, sol.k0 );
    einc = tensor( einc, [ i, q, k, ipol ] );
    hinc = tensor( hinc, [ i, q, k, ipol ] );
    %  Poynting vector in normal direction
    s = 0.5 * dot( nvec,  ...
      sum( w * cross( e - einc, conj( h - hinc ), k ), q ), k );
    %  differential cross section
    ind = indexin( pt, sol.tau );
    dsca( ind{ 1 }, : ) = real( double( s( i, ipol ) ) );    
  end
end

%  impedance and index for propagation direction
Z = arrayfun( @( x ) x.Z( sol.k0 ), obj.layer.mat, 'uniform', 1 );
[ ind1, ind2 ] = deal( obj.dir( :, 3 ) > 0, obj.dir( :, 3 ) <= 0 );
%  normalize differential scattering cross section
dsca( :, ind1 ) = dsca( :, ind1 ) / ( 0.5 / Z( 1 ) );
dsca( :, ind2 ) = dsca( :, ind2 ) / ( 0.5 / Z( end ) );
%  total scattering cross section
csca = sum( dsca, 1 );
