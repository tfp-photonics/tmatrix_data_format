function [ data, ind ] = select2( obj, r, z1, z2, k0, k, path )
%  SELECT2 - Select radii and z-values for given integration path.
%
%  Usage for obj = stratified.isommerfeld :
%    [ data, ind ] = select2( obj, r, z, k0, k, path )
%  Input
%    r      :  radii for reflected Green function
%    z1,z2  :  z-values for reflected Green function
%    k0     :  wavenumber of light in vacuum
%    k      :  wavenumber of light in medium
%    path   :  'semi', 'real' or 'imag'
%  Output
%    data   :  auxiliary structure for evaluation of Green function
%    ind    :  index to selected elements

switch path
  case 'semi'
    data = struct( 'r', r, 'z1', z1, 'z2', z2, 'k0', k0, 'k', k );
    ind = 1 : numel( r );
  case 'real'
    ind = z1 + z2 >= obj.ratio * r;
    data = struct(  ...
      'r', r( ind ), 'z1', z1( ind ), 'z2', z2( ind ), 'k0', k0, 'k', k );
  case 'imag'
    ind = z1 + z2 < obj.ratio * r;
    data = struct(  ...
      'r', r( ind ), 'z1', z1( ind ), 'z2', z2( ind ), 'k0', k0, 'k', k );
end
