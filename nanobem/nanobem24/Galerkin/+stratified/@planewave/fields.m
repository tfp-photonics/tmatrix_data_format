function [ e, h ] = fields( obj, pt, k0 )
%  FIELDS - Electromagnetic fields for planewave excitation.
%
%  Usage for obj = stratified.planewave :
%    [ e, h ] = fields( obj, pt, k0 )
%  Input
%    pt     :  evaluation points
%    k0     :  wavelength of light in vacuum
%  Output
%    e      :  electric field
%    h      :  magnetic field

%  allocate output
[ e, h ] = deal( zeros( numel( pt ), 3, size( obj.pol, 1 ) ) );

%  group points in unique materials
for it = iterpoints( pt )
  
  if it.imat <= obj.layer.n + 1
    %  loop over light propagation directions
    for i1 = 1 : size( obj.dir, 1 )
      %  electric and magnetic fields
      [ e( it.nu, :, i1 ), h( it.nu, :, i1 ) ] = planewave(  ...
        obj.layer, obj.pol( i1, : ), obj.dir( i1, : ), k0, it.pos );
    end  
  end
end
