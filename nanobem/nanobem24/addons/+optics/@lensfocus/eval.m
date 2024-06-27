function field = eval( obj, efield, varargin )
%  EVAL - Evaluate planewave decomposition of focal field.
%
%  Usage for obj = optics.lensfocus :
%    field = eval( obj, efield, PropertyPairs )
%  Input
%    efield     :  incoming electric field before Gaussian sphere
%  PropertyName
%    focus      :  focus position 
%  Output
%    field      :  planewave decomposition of focal fields

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'focus', [ 0, 0, 0 ] );
%  parse input
parse( p, varargin{ : } );

%  wavenumber and refractive index
[ k, n ] = deal( obj.mat.k( obj.k0 ), obj.mat.n( obj.k0 ) );

[ phi, theta ] = deal( obj.phi, obj.theta );
%  unit vectors 
uphi = [ - sin( phi ), cos( phi ), 0 * phi ];
urho = [   cos( phi ), sin( phi ), 0 * phi ];
utheta =  ...
  [ cos( phi ) .* cos( theta ), sin( phi ) .* cos( theta ), - sin( theta ) ];

%  fields after crossing Gaussian reference sphere
efield = bsxfun( @times, uphi,   dot( uphi, efield, 2 ) ) +  ...
         bsxfun( @times, utheta, dot( urho, efield, 2 ) );
%  propagation directions
dir = [ cos( phi ) .* sin( theta ), sin( phi ) .* sin( theta ), cos( theta ) ];

%  prefactor for focused fields
fac = 1i * k * sqrt( 1 / n ) / ( 2 * pi ) * sqrt( cos( theta ) );
%  planewave decomposition of focal fields
efield = bsxfun( @times, efield, fac .* obj.w( : ) .* sin( theta( : ) ) );
field = optics.decompose( obj.k0, efield, dir );

%  manipulate fields after reference sphere
field = field * obj.rot .';
field = field .* exp( - 1i * k * field.dir * p.Results.focus .' );
