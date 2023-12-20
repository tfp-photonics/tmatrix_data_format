function cabs = absorption( obj, sol )
%  ABSORPTION - Absorption cross section.
%
%  Usage for obj = stratified.planewave :
%    cabs = absorption( obj, sol )
%  Input
%    sol    :  solution of BEM equations
%  Output
%    cabs   :  absorption cross section

%  allocate output
cabs = 0;
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
    %  Poynting vector in normal direction
    s = 0.5 * sum( nvec * sum( w * cross( e, conj( h ), k ), q ), i, k );
    cabs = cabs - real( double( s( ipol ) ) );    
  end
end

%  impedance and index for propagation direction
Z = arrayfun( @( x ) x.Z( sol.k0 ), obj.layer.mat, 'uniform', 1 );
[ ind1, ind2 ] = deal( obj.dir( :, 3 ) > 0, obj.dir( :, 3 ) <= 0 );
%  normalize absorption cross section
cabs( :, ind1 ) = cabs( :, ind1 ) / ( 0.5 / Z( 1 ) );
cabs( :, ind2 ) = cabs( :, ind2 ) / ( 0.5 / Z( end ) );
