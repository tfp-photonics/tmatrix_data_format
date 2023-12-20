function [ e, h ] = farfields( obj, dir, k0 )
%  FARFIELD - Electromagnetic far-fields for dipole.
%
%  Usage for obj = stratified.dipole :
%    [ e, h ] = farfields( obj, dir, k0 )
%  Input
%    dir  :  propagation directions for farfields
%    k0   :  wavenumber of light in vacuum
%  Output
%    e    :  electric field
%    h    :  magnetic field

%  material parameters
mat = obj.layer.mat;
data.eps = arrayfun( @( x ) x.eps( k0 ), mat, 'uniform', 1 );
data.mu  = arrayfun( @( x ) x.mu ( k0 ), mat, 'uniform', 1 );
data.Z   = arrayfun( @( x ) x.Z  ( k0 ), mat, 'uniform', 1 );
data.k   = arrayfun( @( x ) x.k  ( k0 ), mat, 'uniform', 1 );
data.k0  = k0;

%  direct contribution
pt = Point( mat, 1, dir );
[ pt( dir( :, 3 ) > 0 ).imat ] = deal( obj.layer.n + 1 );
[ e, h ] = farfields( obj.dip, pt, k0 );
%  avoid directions parallel to z-direction
i1 = abs( dir( :, 3 ) ) == 1;
t = 1e-5;
dir( i1, : ) = dir( i1, : ) *  ...
  [ cos( t ), 0, sin( t ); 0, 1, 0; - sin( t ), 0, cos( t ) ];
%  downgoing and upgoing waves
i1 = dir( :, 3 ) < 0;
i2 = dir( :, 3 ) > 0;

%  loop over source points 
for pt = iterpoints( obj.pt )
  %  source points connected to layer structure ?
  if pt.imat <= obj.layer.n + 1

    %  downgoing waves
    if nnz( i1 )
      [ e1, h1 ] = fields2( obj.layer, data, dir( i1, : ), pt );
      e( i1, :, pt.nu, : ) = e( i1, :, pt.nu, : ) + e1;
      h( i1, :, pt.nu, : ) = h( i1, :, pt.nu, : ) + h1;
    end   
    %  upgoing waves
    if nnz( i2 )
      [ e2, h2 ] = fields2( obj.layer, data, dir( i2, : ), pt );
      e( i2, :, pt.nu, : ) = e( i2, :, pt.nu, : ) + e2;
      h( i2, :, pt.nu, : ) = h( i2, :, pt.nu, : ) + h2;
    end   
  end
end


function [ e, h ] = fields2( layer, data, dir, pt )
%  FIELDS2 - Electromagnetic farfields for source points in single layer.
%
%  Usage :
%    [ e, h ] = fields2( layer, data, dir, pt )
%  Input
%    layer  :  layer structure
%    data   :  auxiliary data
%    pt     :  source points in single layer
%    dir    :  single light propagation direction
%  Output
%    e,h    :  electromagnetic farfield matrices

%  medium index for observation points and interface position
switch sign( dir( 1, 3 ) )
  case 1
    i1 = layer.n + 1;
    z1 = layer.z( end );
  otherwise
    i1 = 1;
    z1 = layer.z( 1 );
end
%  medium index for source points
i2 = pt.imat;
%  wavevector and parallel component
k1 = dir * data.k( i1 );
kpar = sqrt( k1( :, 1 ) .^ 2 + k1( :, 2 ) .^ 2 );
%  number of directions
n = size( dir, 1 );

%  parallel propagation factor and z-values of source points
pos = vertcat( pt.pos );
%  parallel propagation factor and z-values of source points
xpar = k1( :, 1 ) * pos( :, 1 ) .' + k1( :, 2 ) * pos( :, 2 ) .';
z2 = pos( :, 3 );

%  TE and TM basis vectors, Hohenester Eq. (B.5)
te1 = cross( k1, repmat( [ 0, 0, 1 ], n, 1 ), 2 );
te1 = bsxfun( @rdivide, te1, sqrt( dot( te1, te1, 2 ) ) );
tm1 = cross( dir, te1, 2 );

%  coefficients for secondary waves
wave = secondary( layer, data.k0, kpar, i1, i2 );
kz = vertcat( wave.kz );
%  ratio of impedances and wavenumbers
z = data.Z ( i1 ) / data.Z ( i2 );  
r = kz( :, i1 ) ./ kz( :, i2 );
%  allocate output
SL1 = 0;

%  loop over secondary wave propagation directions
for it = 1 : numel( wave( 1 ).te )
  %  propagation direction and distance to interface
  if it == 1 && i2 ~= layer.n + 1
    dir2 = + 1;
    Z = layer.z( i2 ) - z2;
  else
    dir2 = - 1;
    Z = z2 - layer.z( i2 - 1 );
  end
  %  factor from stationary phase approximation
  fac = r .* exp( - 1i * xpar + 1i * kz( :, i2 ) * Z .' ) / ( 4 * pi );
  %  multiply with propagation factor to interface
  fac = fac .* exp( - 1i * k1( :, 3 ) * z1 );
  
  %  wavevector and TE, TM basis vectors
  k2 = [ k1( :, 1 : 2 ), dir2 * kz( :, i2 ) ]; 
  te2 = cross( k2, repmat( [ 0, 0, 1 ], n, 1 ), 2 );
  te2 = bsxfun( @rdivide, te2, sqrt( dot( te2, te2, 2 ) ) );
  tm2 = cross( k2, te2, 2 ) / data.k( i2 );   
  
  %  coefficients for secondary fields
  coeff.te = arrayfun( @( x ) x.te( it ), wave, 'uniform', 1 );
  coeff.tm = arrayfun( @( x ) x.tm( it ), wave, 'uniform', 1 );
  %  dyadic function
  fun = @( x, y, z ) double( tensor( x, [ 1, 2 ] ) *  ...
     tensor( y, [ 1, 3 ] ) * tensor( z, [ 1, 4 ] ), [ 1, 2, 3, 4 ] );
  %  farfield expression of representation formula, Hohenester Eq. (5.37)
  SL1 = SL1 + fun( coeff.te .* te1, fac, te2 ) +  ...
              fun( coeff.tm .* tm1, fac, tm2 ) * z;
end

%  electric dipole farfields
e = data.mu( i2 ) * data.k0 ^ 2 * SL1;
h = double( cross( tensor( dir, [ 1, 2 ] ),  ...
        tensor( e, [ 1, 2, 3, 4 ] ), 2 ), [ 1, 2, 3, 4 ] ) / data.Z( i1 );
