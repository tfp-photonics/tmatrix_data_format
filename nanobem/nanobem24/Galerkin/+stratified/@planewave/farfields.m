function [ e, h, k ] = farfields( obj, k0, key )
%  FARFIELDS - Electromagnetic farfields for planewave excitation.
%
%  Usage for obj = stratified.planewave :
%    [ e, h, k ] = farfields( obj, k0, key )
%  Input
%    k0     :  wavelength of light in vacuum
%    key    :  'refl' or 'trans' for reflected or transmitted wave
%  Output
%    e      :  electric field
%    h      :  magnetic field
%    k      :  wavenumber of reflected or transmitted wave

%  reflected or transmitted wave
if ~exist( 'key', 'var' ),  key = 'refl';  end
%  allocate output
[ e, h, k ] = deal( zeros( size( obj.pol, 1 ), 3 ) );

%  layer structure and small quantity
layer = obj.layer;
eta = 1e-6;
%  wavenumbers
[ k1, z1 ] = deal( layer.mat( 1   ).k( k0 ), layer.z( 1   ) );
[ k2, z2 ] = deal( layer.mat( end ).k( k0 ), layer.z( end ) );
%  loop over light propagation directions
for it = 1 : size( obj.dir, 1 )
  %  light propagation direction and polarization
  [ dir, pol ] = deal( obj.dir( it, : ), obj.pol( it, : ) );
  
  switch sign( dir( 3 ) )
    case 1  %  upgoing wave
      
      %  wavenumbers
      kr = k1 * hypot( dir( 1 ), dir( 2 ) );
      kz1 = stratified.zsqrt( k1 ^ 2 - kr ^  2 );
      kz2 = stratified.zsqrt( k2 ^ 2 - kr ^  2 );
      
      switch key
        case 'refl'
          %  reflected wave for upgoing excitation
          [ er, hr ] = planewave(  ...
            layer, pol, dir, k0, [ 0, 0, z1 - eta ], 'primary', 0 );
          %  multiply with phase factor
          e( it, : ) = er * exp( 1i * kz1 * z1 );
          h( it, : ) = hr * exp( 1i * kz1 * z1 );
          k( it, : )= [ kr * dir( 1 ), kr * dir( 2 ), - kz1 ];
        case 'trans'
          %  transmitted wave for upgoing excitation
          [ et, ht ] = planewave(  ...
            layer, pol, dir, k0, [ 0, 0, z2 + eta ], 'primary', 0 );
          %  multiply with phase factor         
          if kr < k2
            e( it, : ) = et * exp( - 1i * kz2 * z2 );
            h( it, : ) = ht * exp( - 1i * kz2 * z2 );
            k( it, : )= [ kr * dir( 1 ), kr * dir( 2 ), kz2 ];
          end
      end
      
    otherwise  %  downgoing wave
      
      %  wavenumbers
      kr = k2 * hypot( dir( 1 ), dir( 2 ) );
      kz1 = stratified.zsqrt( k1 ^ 2 - kr ^  2 );
      kz2 = stratified.zsqrt( k2 ^ 2 - kr ^  2 );
      
      switch key
        case 'refl'
          %  reflected wave for downgoing excitation
          [ er, hr ] = planewave(  ...
            layer, pol, dir, k0, [ 0, 0, z2 + eta ], 'primary', 0 );
          %  multiply with phase factor
          e( it, : ) = er * exp( - 1i * kz2 * z2 );
          h( it, : ) = hr * exp( - 1i * kz2 * z2 );
          k( it, : )= [ kr * dir( 1 ), kr * dir( 2 ), kz2 ];
        case 'trans'
          %  transmitted wave for downgoing excitation
          [ et, ht ] = planewave(  ...
            layer, pol, dir, k0, [ 0, 0, z1 - eta ], 'primary', 0 );
          %  multiply with phase factor         
          if kr < k1
            e( it, : ) = et * exp( 1i * kz1 * z1 );
            h( it, : ) = ht * exp( 1i * kz1 * z1 );
            k( it, : )= [ kr * dir( 1 ), kr * dir( 2 ), - kz1 ];
          end
      end      
  end
end

