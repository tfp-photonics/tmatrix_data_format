function [ data, ind ] = select( obj, r, z, k0, k, path )
%  SELECT - Select radii and z-values for given integration path.
%
%  Usage for obj = stratified.isommerfeld :
%    [ data, ind ] = select( obj, r, z, k0, k, path )
%  Input
%    r,z    :  radii and z-values for reflected Green function
%    k0     :  wavenumber of light in vacuum
%    k      :  wavenumber of light in medium
%    path   :  'semi', 'real' or 'imag'
%  Output
%    data   :  auxiliary structure for evaluation of Green function
%    ind    :  index to selected elements

switch class( z )
  case 'struct'
    [ data, ind ] = select2( obj, r, z.z1, z.z2, k0, k, path );
  otherwise
    [ data, ind ] = select1( obj, r, z, k0, k, path );
end
