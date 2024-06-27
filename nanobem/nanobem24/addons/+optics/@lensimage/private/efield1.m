function e = efield1( obj, far, x, y, varargin )
%  EFIELD1 - Electric field on image side, Eq. (3.10).
%
%  Usage for obj = optics.lensimage :
%    e = efield1( obj, far, x, y, PropertyPairs )
%  Input
%    far    :  electric far-field along directions of obj.dir 
%    x,y    :  image coordinates
%  PropertyName
%    focus  :  focus position of imaging lens
%    ntab   :  number of tabulation points for Bessel function
%  Output
%     e     :  electric image fields in image plane

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'focus', [ 0, 0, 0 ] );
addParameter( p, 'ntab', 2001 );
addParameter( p, 'waitbar', 0 );
%  parse input
parse( p, varargin{ : } );

%  refractive indices and wavenumbers
n1 = obj.mat1.n( obj.k0 );  k1 = obj.k0 * n1;
n2 = obj.mat2.n( obj.k0 );  k2 = obj.k0 * n2;
%  manipulate fields before reference sphere
far = far .* exp( 1i * k1 * obj.dir * p.Results.focus .' );
far = far * obj.rot;

[ up1, ut1, ur1 ] = optics.cart2unit( obj.dir * obj.rot );
%  transformation for change of coordinate systems
far = optics.unitproject( far, ut1, ur1 ) +  ...
      optics.unitproject( far, up1, up1 );
%  apply manipulation function in backfocal plane
if ~isempty( obj.backfocal ), far = obj.backfocal( far );  end   

%  polar coordinates on image side
[ xx, yy ] = ndgrid( x, y );
[ u2, r2 ] = cart2pol( xx( : ), yy( : ) );

%  number of angular coordinates
m1 = numel( obj.theta );
m2 = numel( obj.phi );
%  Fourier transform for azimuthal angle
mtab = - fix( 0.5 * m2 ) : fix( 0.5 * m2 );
far = fftshift( fft( reshape( far, m2, m1, 3, [] ), [], 1 ), 1 ) / m2;
%  argument for Bessel functions
z = k1 * r2 * sin( obj.theta );
ztab = linspace( min( z( : ) ), max( z( : ) ), p.Results.ntab );

%  prefactor for Richards-Wolf integral
[ sint, cost ] = deal( sin( obj.theta ), cos( obj.theta ) );
fac = 1i * k2 * ( n1 / n2 ) .^ 1.5 / ( 2 * pi ) * obj.w .* sqrt( cost ) .* sint;
%  allocate output
e = 0;

if p.Results.waitbar,  multiWaitbar( 'lensimage', 0, 'Color', 'g' );  end
%  loop over angular degrees
for m = 0 : obj.mcut
  %  Bessel function 
  j = interp1( ztab, besselj( m, ztab ), z );
  %  perform integration using Eq. (3.21)
  rho = 2 * pi * 1i ^ m * j .* ( exp( + 1i * m * u2 ) * fac );
  e = e + rho * reshape( far( + m == mtab, :, :, : ), m1, [] );
  %  deal with negative orders
  if m ~= 0 
    rho = 2 * pi * 1i ^ m * j .* ( exp( - 1i * m * u2 ) * fac );
    e = e + rho * reshape( far( - m == mtab, :, :, : ), m1, [] );
  end   
  
  if p.Results.waitbar && ~mod( m, 10 )
    multiWaitbar( 'lensimage', m / obj.mcut );  
  end
end
%  close waitbar
if p.Results.waitbar,  multiWaitbar( 'CloseAll' );  end

%  reshape output
e = reshape( e, numel( x ), numel( y ), 3, [] );
