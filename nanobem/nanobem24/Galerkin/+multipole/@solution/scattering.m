function sca = scattering( obj, varargin )
%  SCATTERING - Scattered power, Hohenester Eq. (E.29).
%
%  Usage for obj = multipole.solution :
%    sca = scattering( obj, PropertyPairs )
%  PropertyName
%    pinfty   :  sphere segment for far-fields
%  Output
%    sca      :  scattered power

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'pinfty', [] );
addParameter( p, 'shift', [] );
%  parse input
parse( p, varargin{ : } );

if isempty( p.Results.pinfty )
  %  wavenumber and impedance
  [ k, Z ] = deal( obj.mat.k( obj.k0 ), obj.mat.Z( obj.k0 ) );  
  %  scattered power, Hohenester Eq. (E.29)
  sca = 0.5 * Z / k ^ 2 * sum( abs( obj.a ) .^ 2 + abs( obj.b ) .^ 2 ); 
else
  %  electromagnetic far-fields
  pinfty = p.Results.pinfty;
  [ e, h ] = farfields( obj, pinfty.pos, 'shift', p.Results.shift );
  %  scattering cross section
  siz = size( e );
  s = 0.5 * real( cross( e, conj( h ), 2 ) );
  sca = sum( reshape( pinfty.area .* pinfty.pos, [], 1 ) .*  ...
                          reshape( s, siz( 1 ) * 3, [] ), 1 );
  %  reshape output
  sca = reshape( sca, [ siz( 3 : end ), 1 ] );
end
