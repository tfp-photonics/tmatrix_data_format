function [ f, n, data ] = optforce( obj, varargin )
%  OPTFORCE - Optical force and torque.
%
%  Usage for obj = multipole.solution :
%    [  f, n, data ] = optforce( obj, data )
%  Input
%    data   :  auxiliary data for computation of force
%  Output
%    f      :  optical force  in pN
%    n      :  optical torque in pN × nm
%    data   :  auxiliary data

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'data', [] );
%  parse input
parse( p, varargin{ : } );

if isempty( p.Results.data )
  data = init( obj );
else
  data = p.Results.data;
end

%  Mie coefficients, account for different definitions
tab = obj.tab;
fac = ( - 1i ) .^ ( tab.l + 1 );
[ a, ai ] = deal( fac .* obj.b, fac .* obj.bi );
[ b, bi ] = deal( fac .* obj.a, fac .* obj.ai );
%  allocate output
[ f, n ] = deal( zeros( size( a, 2 ), 3 ) );

%  loop over Cartesian coordinates
for k = 1 : 3
  %  evaluation functions
  fun1 = fun( data{ k }.I1 );  
  fun2 = fun( data{ k }.I2 );
  fun3 = fun( data{ k }.I3 );
  fun4 = fun( data{ k }.I4 );
  %  optical force
  f( :, k ) = fun1( a, a ) + fun1( b, b ) + fun1( ai, a ) + fun1( bi, b ) -  ...
              fun2( b, a ) + fun2( a, b ) - fun2( bi, a ) + fun2( ai, b );
  % optical torque
  n( :, k ) = fun3( a, a ) + fun3( b, b ) + fun3( ai, a ) + fun3( bi, b ) -  ...
              fun4( a, b ) + fun4( b, a );
end

%  wavenumber and permeability of embedding medium
mat = obj.mat;
[ k, mu ] = deal( mat.k( obj.k0 ), mat.mu( obj.k0 ) );
%  conversion factor force in pN, use vacuum permittivity
fac = 1e12 * 8.854e-12;
%  force and torque
f = - 0.5 * mu * fac / k ^ 2 * real( f );
n = - 0.5 * mu * fac / k ^ 3 * real( n );



function fun = fun( mat )
%  FUN - Evalulation function

%  nonzero elements of matrix
[ i1, i2, val ] = find( mat );
%  anonymous function
fun = @( a, b ) sum( conj( a( i1, : ) ) .* val .* b( i2, : ), 1 );


function data = init( obj )
%  INIT - Precompute coefficients for computation of optical forces.

%  quadrature points and weights
quad = quadsph( obj, 'vector' );
%  grid for spherical harmonics
tab = obj.tab;
n = numel( tab.l );
[ i1, i2 ] = ndgrid( 1 : n );
%  angular degree and orders changes by one
ind = abs( tab.l( i1 ) - tab.l( i2 ) ) <= 1 &  ...
      abs( tab.m( i1 ) - tab.m( i2 ) ) <= 1;
[ i1, i2 ] = deal( i1( ind ), i2( ind ) );

%  vector spherical harmonics and spherical harmonics
[ xm, xe, ~, y ] = vsh( tab.l, tab.m, quad.t, quad.u );
xe = 1i * xe;
y = y .* sqrt( tab.l .* ( tab.l + 1 ) );

%  unit vector
[ u, t ] = deal( quad.u, quad.t );
pos = [ cos( u ) .* sin( t ), sin( u ) .* sin( t ), cos( t ) ];
%  allocate output
data = cell( 1, 3 );

%  loop over Cartesian coordinates
for k = 1 : 3
  %  coefficients for force evaluation
  z = repmat( pos( :, k ) .* quad.w, 3, 1 );
  data{ k }.I1 = dot( xm( i1, : ), xm( i2, : ) .* z .', 2 );
  data{ k }.I2 = dot( xe( i1, : ), xm( i2, : ) .* z .', 2 );
  %  coefficients for torque evaluation
  data{ k }.I3 = dot( y( i1, : ), xm( i2, :, k ) .* quad.w .', 2 );
  data{ k }.I4 = dot( y( i1, : ), xe( i2, :, k ) .* quad.w .', 2 );
  %  assemble and compress arrays
  for name = [ "I1", "I2", "I3", "I4" ]
    z = data{ k }.( name );
    z( abs( z ) < 1e-10 ) = 0;
    data{ k }.( name ) = accumarray( { i1, i2 }, z, [ n, n ], [], [], true );
  end
end
