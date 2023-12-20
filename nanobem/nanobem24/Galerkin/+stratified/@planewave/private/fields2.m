function [ e, h ] = fields2( obj, pt, k0 )
%  FIELDS2 - Electromagnetic fields for planewave excitation.
%
%  Usage for obj = stratified.planewave :
%    [ e, h ] = fields2( obj, pt, k0 )
%  Input
%    pt     :  integration points
%    k0     :  wavelength of light in vacuum
%  Output
%    e      :  electric field
%    h      :  magnetic field

%  positions 
pos = eval( pt );
%  allocate output
e = zeros( [ size( pos ), size( obj.pol, 1 ) ] );
h = zeros( [ size( pos ), size( obj.pol, 1 ) ] );

if pt.inout( 2 ) <= obj.layer.n + 1
  %  loop over light propagation directions
  for it = 1 : size( obj.dir, 1 )
    %  electric and magnetic fields
    [ e1, h1 ] = planewave( obj.layer,  ...
      obj.pol( it, : ), obj.dir( it, : ), k0, reshape( pos, [], 3 ) );
    %  accumulate electromagnetic fields
    e( :, :, :, it ) = reshape( e1, size( pos ) );
    h( :, :, :, it ) = reshape( h1, size( pos ) );
  end
end
