function wave = secondary( obj, k0, kpar, i1, i2, varargin )
%  SECONDARY - Coefficients for secondary waves.
%
%  Usage for obj = stratified.layerstructure :
%    wave = secondary( obj, k0, kpar, i1, i2 )
%  Input
%    k0     :  wavenumber of light in vacuum
%    kpar   :  parallel momenta
%    i1     :  medium index for observation point
%    i2     :  medium index for source point
%  Output
%    wave   :  coefficients for secondary waves

n = obj.n;
%  transfer matrices
obj = eval( obj, k0 );
m = arrayfun( @( it ) transfer(  ...
  obj, k0, kpar, it, varargin{ : } ), 1 : n, 'uniform', 1 );
%  wavenumbers
fun = @( kpar ) stratified.zsqrt( obj.data.eps .* obj.data.mu * k0 ^ 2 - kpar ^ 2 );
kz = arrayfun( fun, kpar, 'uniform', 0 );
wave = struct( 'kz', kz );

%  media indices for up- and downgoing waves
%    only one wave in media outside of layer structure
ind = repmat( 1 : n + 1, 2, 1 );
ind = ind( 2 : end - 1 );
%  index for connected primary and secondary fields
ind1 = accumarray( [ 1; 2 * n + 2 ], false, [ 2 * n + 2, 1 ], @( x ) x, true );
ind2 = accumarray( [ 2; 2 * n + 1 ], false, [ 2 * n + 2, 1 ], @( x ) x, true );
%  initialize block matrices
[ lhs, rhs ] = deal( zeros( 2 * n, 2 * n + 2 ) );

%  loop over parallel momenta
for ipar = 1 : numel( kpar )

  %  propagation constants
  fac = exp( 1i * [ 0, wave( ipar ).kz( 2 : n ) .* obj.d( 1 : n - 1 ), 0 ] );
  %  set block matrices to zero
  lhs = 0 * lhs;
  rhs = 0 * rhs;

  for name = [ "te", "tm" ]
    %  fill matrices
    for it = 1 : n
      %  transfer matrix and indices for sub-matrices
      M = m( it ).( name )( :, :, ipar );
      k1 = 2 * it - 2 + ( 1 : 2 );
      k2 = 2 * it     + ( 1 : 2 );
      %  fill sub-matrices
      lhs( k1, k1 ) = [ fac( it ), 0; 0, 1 ];
      lhs( k1, k2 ) = - M * [ 1, 0; 0, fac( it + 1 ) ];
      rhs( k1, k1 ) = [ 1, 0; 0, 0 ];
      rhs( k1, k2 ) = - M * [ 0, 0; 0, 1 ];
    end
     %  secondary fields
    x = - lhs( :, ind1 ) \ rhs( :, ind2 ); 
    wave( ipar ).( name ) = x( ind == i1, ind == i2 );
  end
end

