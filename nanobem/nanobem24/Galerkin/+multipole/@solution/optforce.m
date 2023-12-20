function [ f, data ] = optforce( obj, varargin )
%  OPTFORCE - Optical force.
%    See Guiterrez-Cuevas et al., PRA 97, 053448 (2018).
%
%  Usage for obj = multipole.solution :
%    [ f, data ] = optforce( obj, data )
%  Input
%    data   :  auxiliary data for computation of force
%  Output
%    f      :  optical force in pN
%    data   :  auxiliary data for reuse

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

%  wavenumber in embedding medium
k2 = obj.mat.k( obj.k0 );
%  Mie coefficients, account for different definitions 
a  = @( i ) obj.b( i, : );  
b  = @( i ) obj.a( i, : );
ai = @( i ) obj.bi( i, : );
bi = @( i ) obj.ai( i, : );

%  auxiliary multiplication function
fun = @( a, b ) conj( a ) .* b;

%  first term in Eq. (29a)
[ i1, i2, val ] = find( data.fx1 );
fx1 = 2 * fun( a( i1 ), a( i2 ) ) + fun( ai( i1 ), a( i2 ) ) + fun( a( i1 ), ai( i2 ) ) +  ...
      2 * fun( b( i1 ), b( i2 ) ) + fun( bi( i1 ), b( i2 ) ) + fun( b( i1 ), bi( i2 ) );
fx1 = sum( bsxfun( @times, fx1, val ), 1 );
%  second term in Eq. (29a)
[ i1, i2, val ] = find( data.fx2 );
fx2 = 2 * fun( a( i1 ), a( i2 ) ) + fun( ai( i1 ), a( i2 ) ) + fun( a( i1 ), ai( i2 ) ) +  ...
      2 * fun( b( i1 ), b( i2 ) ) + fun( bi( i1 ), b( i2 ) ) + fun( b( i1 ), bi( i2 ) );
fx2 = sum( bsxfun( @times, fx2, val ), 1 );
%  third term in Eq. (29a)
[ i1, i2, val ] = find( data.fx3 );
fx3 = 2 * fun( a( i1 ), b( i2 ) ) + fun( ai( i1 ), b( i2 ) ) + fun( a( i1 ), bi( i2 ) ) -  ...
      2 * fun( b( i1 ), a( i2 ) ) - fun( bi( i1 ), a( i2 ) ) - fun( b( i1 ), ai( i2 ) );
fx3 = sum( bsxfun( @times, fx3, val ), 1 );    

%  total force in x and y directions
fx = 0.25i / k2 ^ 2 * ( fx1 + fx2 + fx3 );

%  first term in Eq. (29b)
[ i1, i2, val ] = find( data.fz1 );
fz1 = 2 * fun( a( i1 ), a( i2 ) ) + fun( ai( i1 ), a( i2 ) ) + fun( a( i1 ), ai( i2 ) ) +  ...
      2 * fun( b( i1 ), b( i2 ) ) + fun( bi( i1 ), b( i2 ) ) + fun( b( i1 ), bi( i2 ) );
fz1 = sum( bsxfun( @times, fz1, val ), 1 );
%  second term in Eq. (29b)
[ i1, i2, val ] = find( data.fz2 );
fz2 = 2 * fun( b( i1 ), a( i2 ) ) + fun( bi( i1 ), a( i2 ) ) + fun( b( i1 ), ai( i2 ) );
fz2 = sum( bsxfun( @times, fz2, val ), 1 );

%  total force in z-direction
fz = - 0.5 / k2 ^ 2 * imag( fz1 + fz2 );
%  assemble force
f = [ real( fx( : ) ), imag( fx( : ) ), fz( : ) ];
%  conversion factor force in pN, use vacuum permittivity
fac = 1e12 * 8.854e-12;
f = fac * f;


function data = init( obj )
%  INIT - Precompute coefficients for computation of optical forces.
%    See Guiterrez-Cuevas et al., PRA 97, 053448 (2018), Eq. (29).

%  table of angular degrees and orders
tab = [ obj.tab.l( : ), obj.tab.m( : ) ];
%  allocate output
[ fx1, fx2, fx3, fz1, fz2 ] = deal( zeros( size( tab, 1 ) ) );

%  loop over angular degrees and orders
for l = 1 : max( obj.tab.l )
for m = - unique( tab( tab( :, 1 ) == l, : ) ) .'
  %  first term in Eq. (29a)
  %    prefactors may differ because of different function definitions
  [ ~, i1 ] = ismember( [ l + 1, m + 1 ], tab, 'rows' );
  [ ~, i2 ] = ismember( [ l,     m     ], tab, 'rows' );
  if i1 && i2
    r1 = ( l + m + 2 ) * ( l + m + 1 ) * l * ( l + 2 );
    r2 = ( 2 * l + 1 ) * ( 2 * l + 3 );
    fx1( i1, i2 ) = 1 / ( l + 1 ) * sqrt( r1 / r2 );
  end
  %  second term in Eq. (29a)
  [ ~, i1 ] = ismember( [ l,     m     ], tab, 'rows' );
  [ ~, i2 ] = ismember( [ l + 1, m - 1 ], tab, 'rows' );
  if i1 && i2
    r1 = ( l - m + 2 ) * ( l - m + 1 ) * l * ( l + 2 );
    r2 = ( 2 * l + 1 ) * ( 2 * l + 3 );
    fx2( i1, i2 ) = 1 / ( l + 1 ) * sqrt( r1 / r2 );
  end  
  %  third term in Eq. (29a)
  [ ~, i1 ] = ismember( [ l, m + 1 ], tab, 'rows' );
  [ ~, i2 ] = ismember( [ l, m     ], tab, 'rows' ); 
  if i1 && i2
    r1 = ( l + m + 1 ) * ( l - m );
    r2 = l * ( l + 1 );
    fx3( i1, i2 ) = - sqrt( r1 ) ./ r2;
  end    
  
  %  first term in Eq. (29b)
  [ ~, i1 ] = ismember( [ l,     m ], tab, 'rows' );
  [ ~, i2 ] = ismember( [ l + 1, m ], tab, 'rows' );
  if i1 && i2 
    r1 = ( l - m + 1 ) * ( l + m + 1 ) * l * ( l + 2 );
    r2 = ( 2 * l + 1 ) * ( 2 * l + 3 );
    fz1( i1, i2 ) = 1 / ( l + 1 ) * sqrt( r1 / r2 );
  end
  %  second term in Eq. (29b)
  [ ~, i1 ] = ismember( [ l, m ], tab, 'rows' );
  fz2( i1, i1 ) = m / ( l * ( l + 1 ) );
end
end

%  set output
data = struct(  ...
  'fx1', fx1, 'fx2', fx2, 'fx3', fx3, 'fz1', fz1, 'fz2', fz2 );
