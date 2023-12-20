function q = qinc( obj, tau, k0, varargin )
%  QINC - Inhomogeneities for multipole excitation.
%
%  Usage for obj = multipole.tsolver :
%    q = qinc( obj, tau, k0, PropertyPairs )
%  Input
%    tau      :  discretized particle boundary
%    k0       :  wavenumber of light in vacuum
%  Output
%    q        :  multipole inhomogeneities for BEM solver

%  wavenumber of embedding medium
mat = obj.embedding;
[ k1, Z1 ] = deal( mat.k( k0 ), mat.Z( k0 ) );

%  dummy indices for tensor class
[ i, q, a, k, n ] = deal( 1, 2, 3, 4, 5 );
%  allocate multipole coefficients
n1 = numel( obj.tab.l );
[ M, N ] = deal( zeros( ndof( tau ), n1 ) );

%  loop over boundary elements
for pt = quadboundary( tau, obj.rules )
  
  if any( pt.tau( 1 ).inout == obj.imat )
    %  quadrature points, weights, and shape functions
    [ pos, w, f ] = eval( pt );
    w = tensor( w, [ i, q ] );
    f = tensor( f, [ i, q, a, k ] );
    
    %  transverse vector functions
    [ M1, N1 ] = transvec( obj, reshape( pos, [], 3 ), k1 );
    %  convert to tensor class
    M1 = tensor( reshape( M1, [ n1, size( pos ) ] ), [ n, i, q, k ] );
    N1 = tensor( reshape( N1, [ n1, size( pos ) ] ), [ n, i, q, k ] );        
    %  perform integrations
    M1 = double( sum( w * dot( M1, f, k ), q ), [ i, a, n ] );
    N1 = double( sum( w * dot( N1, f, k ), q ), [ i, a, n ] );
    %  reshape arrays
    M1 = reshape( M1, [], size( M1, 3 ) );
    N1 = reshape( N1, [], size( N1, 3 ) );
    
    %  add to matrices
    nu = vertcat( pt.tau.nu );
    %  loop over global degrees of freedom  
    for it = 1 : numel( nu )
      M( nu( it ), : ) = M( nu( it ), : ) + M1( it, : );
      N( nu( it ), : ) = N( nu( it ), : ) + N1( it, : );
    end        
  end
end

%  multipoles, see Hohenester "Nano and Quantum Optics", Eq. (E.4)
e = [ N,   M ] * Z1;
h = [ M, - N ]; 
%  set output
q = struct( 'e', e, 'h', h, 'tau', tau, 'k0', k0 );

