function e = efield2( obj, field, x, y, varargin )
%  EFIELD2 - Electric field on image side for planewave decomposition.
%
%  Usage for obj = optics.lensimage :
%    e = efield2( obj, field, x, y, PropertyPairs )
%  Input
%    field  :  planewave decomposition
%    x,y    :  image coordinates
%  PropertyName
%    focus  :  focus position of imaging lens
%  Output
%     e     :  electric image fields in image plane

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'focus', [ 0, 0, 0 ] );
%  parse input
parse( p, varargin{ : } );

%  refractive indices
n1 = obj.mat1.n( obj.k0 );  k1 = obj.mat1.k( obj.k0 );
n2 = obj.mat2.n( obj.k0 );  
%  manipulate fields before reference sphere
field = field .* exp( + 1i * k1 * field.dir * p.Results.focus .' );
field = field * obj.rot;

%  normalization function
norm = @( x ) bsxfun( @rdivide, x, sqrt( dot( x, x, 2 ) ) );
%  unit vectors 
ez  = repmat( [ 0, 0, 1 ], field.n, 1 );
te  = norm( cross( field.dir, ez, 2 ) );  %  up1
if any( isnan( te( : ) ) )
  [ te( isnan( te( :, 1 ) ), : ) ] = deal( [ 1, 0, 0 ] );
end
tm1 = norm( cross( field.dir, te, 2 ) );  %  ut1
tm2 = norm( cross( ez, te, 2 ) );         %  ur1

%  prefactor for electric field
theta = acos( field.dir( :, 3 ) );
fac = sqrt( n1 / n2 ) * sqrt( cos( theta ) );
%  propagate electric field through lens
%    Khadir et al., JOSA 36, 478 (2019), Eqs. (2,3)
efield = bsxfun( @times, te,  fac .* dot( te,  field.efield, 2 ) ) +  ...
         bsxfun( @times, tm2, fac .* dot( tm1, field.efield, 2 ) );
%  apply manipulation function in backfocal plane
if ~isempty( obj.backfocal ), efield = obj.backfocal( efield );  end
       
%  coordinates on image side 
[ xx, yy ] = ndgrid( x, y );
pos = [ xx( : ), yy( : ), 0 * xx( : ) ];
%  wavenumber in medium and allocate output
k1 = obj.k0 * n1;  
e = 0;
%  loop over wavevectors inside acceptance cone of NA
for i1 = find( theta <= max( obj.theta ) ) .'
  e = e + exp( 1i * k1 * pos * field.dir( i1, : ) .' ) * efield( i1, : );
end
%  reshape output
e = reshape( e, [ size( xx ), 3 ] );
