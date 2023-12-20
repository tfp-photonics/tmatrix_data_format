function [ e, h ] = fields2( obj, data, dir, pt )
%  FIELDS2 - Electromagnetic farfields for source points in single layer.
%
%  Usage for obj = stratified.solution :
%    [ e, h ] = fields2( obj, data, dir, pt )
%  Input
%    data   :  auxiliary data
%    pt     :  quadrature points for source points in single layer
%    dir    :  light propagation directions
%  Output
%    e,h    :  electromagnetic farfields

%  medium index for observation points and interface position
switch sign( dir( 1, 3 ) )
  case 1
    i1 = obj.layer.n + 1;
    z1 = obj.layer.z( end );
  otherwise
    i1 = 1;
    z1 = obj.layer.z( 1 );
end
%  medium index for source points
i2 = pt.tau( 1 ).inout( 2 );
%  wavevector and parallel component
k1 = dir * data.k( i1 );
kpar = sqrt( k1( :, 1 ) .^ 2 + k1( :, 2 ) .^ 2 );
%  number of directions
n = size( dir, 1 );

%  dummy indices for tensor class
[ i, nu, q, a, k ] = deal( 1, 2, 3, 4, 5 );
%  parallel propagation factor and z-values of source points
xpar = tensor( k1( :, 1 ), nu ) * tensor( data.pos( :, :, 1 ), [ i, q ] ) +  ...
       tensor( k1( :, 2 ), nu ) * tensor( data.pos( :, :, 2 ), [ i, q ] );
z2 = tensor( data.pos( :, :, 3 ), [ i, q ] );
%  convert shape elements to tensor class
f = tensor( data.f, [ i, q, a, k ] ) * tensor( data.w, [ i, q ] );

%  TE and TM basis vectors, Hohenester Eq. (B.5)
te1 = cross( k1, repmat( [ 0, 0, 1 ], n, 1 ), 2 );
te1 = bsxfun( @rdivide, te1, sqrt( dot( te1, te1, 2 ) ) );
tm1 = cross( dir, te1, 2 );
%  convert to tensor class
te1 = tensor( te1, [ nu, k ] );
tm1 = tensor( tm1, [ nu, k ] );

%  coefficients for secondary waves 
wave = secondary( obj.layer, obj.k0, kpar, i1, i2 );
kz = vertcat( wave.kz );
%  ratio of impedances and wavenumbers
z = data.Z( i1 ) / data.Z( i2 );  
r = tensor( kz( :, i1 ) ./ kz( :, i2 ), nu );
%  allocate output
[ SL1, DL1 ] = deal( 0 );

%  loop over secondary wave propagation directions
for it = 1 : numel( wave( 1 ).te )
  %  propagation direction and distance to interface
  if it == 1 && i2 ~= obj.layer.n + 1
    dir2 = + 1;
    Z = obj.layer.z( i2 ) - z2;
  else
    dir2 = - 1;
    Z = z2 - obj.layer.z( i2 - 1 );
  end
  %  factor from stationary phase approximation
  kz2 = tensor( kz( :, i2 ), nu );
  fac = r * exp( - 1i * xpar + 1i * kz2 * Z ) / ( 4 * pi );  
  %  multiply with propagation factor to interface
  kz1 = tensor( k1( :, 3 ), nu ); 
  fac = fac * exp( - 1i * kz1 * z1 );
  
  %  wavevector and TE, TM basis vectors 
  k2 = [ k1( :, 1 : 2 ), dir2 * kz( :, i2 ) ];
  te2 = cross( k2, repmat( [ 0, 0, 1 ], n, 1 ), 2 );
  te2 = bsxfun( @rdivide, te2, sqrt( dot( te2, te2, 2 ) ) );
  tm2 = cross( k2, te2, 2 ) / data.k( i2 );
  %  convert to tensor class
  te2 = tensor( te2, [ nu, k ] );
  tm2 = tensor( tm2, [ nu, k ] );
    
  %  assemble function
  fun = @( x ) reshape( double( x, [ nu, k, i, a ] ), n, 3, [] );
  %  coefficients for secondary fields
  coeff.te = tensor( arrayfun( @( x ) x.te( it ), wave, 'uniform', 1 ), nu );
  coeff.tm = tensor( arrayfun( @( x ) x.tm( it ), wave, 'uniform', 1 ), nu );
  %  farfield expression of representation formula, Hohenester Eq. (5.37)
  fte2 = sum( fac * sum( te2 * f, k ), q );
  ftm2 = sum( fac * sum( tm2 * f, k ), q );
  SL1 = SL1 + fun( coeff.te * te1 * fte2 + coeff.tm * tm1 * ftm2 * z );
  DL1 = DL1 + fun( coeff.te * te1 * ftm2 - coeff.tm * tm1 * fte2 * z );
end

%  multiplication function
fun = @( x, y ) reshape( reshape( x, [], size( x, 3 ) ) * y, n, 3, [] );
%  put together electomagnetic fields
e = 1i * data.mu ( i2 ) * obj.k0 * fun( SL1, data.uh ) -   ...
               1i * data.k( i2 ) * fun( DL1, data.ue );
h = double( cross( tensor( dir, [ nu, k ] ),  ...
           tensor( e, [ nu, k, i ] ), k ), [ nu, k, i ] ) / data.Z( i1 );
