function t = eval( obj, sol, varargin )
%  EVAL - Multipole coefficients from BEM solution.
%
%  Usage for obj = multipole.tsolver :
%    t = eval( sol, PropertyPairs )
%  Input
%    sol  :  BEM solution
%  PropertyName
%    mat  :  compute T-matrix of multipole coefficients
%  Output
%    t    :  T-matrix or multipole coefficients

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'mat', size( sol.e, 2 ) == 2 * numel( obj.tab.l ) );
%  parse input
parse( p, varargin{ : } );

%  discretized particle and wavenumber of light in vacuum
[ tau, k0 ] = deal( sol.tau, sol.k0 );
%  wavenumber and impedance of embedding medium
mat = obj.embedding;
[ k1, Z1 ] = deal( mat.k( k0 ), mat.Z( k0 ) );
%  spherical degrees
ltab = obj.tab.l;
%  spherical degrees and number of solution vectors
[ n, m ] = deal( numel( ltab ), size( sol.e, 2 ) );

%  dummy indices for tensor class
[ i, j, q, k ] = deal( 1, 2, 3, 4 );
%  allocate multipole coefficients
[ a, b ] = deal( zeros( numel( ltab ), size( sol.e, 2 ) ) );

%  loop over boundary elements
for pt = quadboundary( tau, obj.rules )
  
  if any( pt.tau( 1 ).inout == obj.imat )
    %  quadrature points and weights
    [ pos, w ] = eval( pt );
    w = tensor( w, [ i, q ] );
    %  normal vector
    nvec = tensor( vertcat( pt.tau.nvec ), [ i, k ] );
    
    %  transverse vector functions
    [ M, N ] = transvec( obj, reshape( pos, [], 3 ), k1, 'conj', 1 );
    %  interpolate tangential fields
    [ e, h ] = interp( sol, pt );
    ue = double( w * cross( nvec, tensor( e, [ i, q, k, j ] ), k ), [ i, q, k, j ] );
    uh = double( w * cross( nvec, tensor( h, [ i, q, k, j ] ), k ), [ i, q, k, j ] );
    
    %  perform integrations
    Me = reshape( M, n, [] ) * reshape( ue, [], m );
    Mh = reshape( M, n, [] ) * reshape( uh, [], m );
    Ne = reshape( N, n, [] ) * reshape( ue, [], m );
    Nh = reshape( N, n, [] ) * reshape( uh, [], m );
    %  update multipole coefficients
    a = a + k1 ^ 2 * ( Me / Z1 - Nh );
    b = b - k1 ^ 2 * ( Ne / Z1 + Mh );
  end
end

%  set output
switch p.Results.mat
  case 1
    n = numel( ltab );
    t = multipole.tmatrix( obj, k0 );
    t.aa = a( :, 1 : n );  t.ab = a( :, n + 1 : end );
    t.ba = b( :, 1 : n );  t.bb = b( :, n + 1 : end );
  otherwise
    t = multipole( obj.tab, mat, k0, a, b );
end
