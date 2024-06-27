function [ polx, poly ] = polarization( obj, u )
%  POLARIZATION - Linear, radial or azimuthal polarization vector(s).
%
%  Usage for obj = laguerregauss :
%    [ polx, poly ] = polarization( obj, u )
%  Input
%    u      :  azimuthal angle
%  Output
%    polx   :  polarization vector along x
%    poly   :  polarization vector along y

if isnumeric( obj.pol )
  %  assert that incoming polarization is perpendicular to z-direction
  assert( obj.pol( 3 ) == 0 );
  [ polx, poly ] = deal( obj.pol( 1 ), obj.pol( 2 ) );
else
  switch obj.pol( 1 : 3 )
    case 'rad'
      [ polx, poly ] = deal(   cos( u ), sin( u ) );
    case 'azi'
      [ polx, poly ] = deal( - sin( u ), cos( u ) );
  end
end
