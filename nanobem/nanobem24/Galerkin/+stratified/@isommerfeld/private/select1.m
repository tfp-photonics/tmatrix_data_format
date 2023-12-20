function [ data, ind ] = select1( obj, r, z, k0, k, path )
%  SELECT1 - Select radii and z-values for given integration path.
%
%  Usage for obj = stratified.isommerfeld :
%    [ data, ind ] = select1( obj, r, z, k0, k, path )
%  Input
%    r,z    :  radii and z-values for reflected Green function
%    k0     :  wavenumber of light in vacuum
%    k      :  wavenumber of light in medium
%    path   :  'semi', 'real' or 'imag'
%  Output
%    data   :  auxiliary structure for evaluation of Green function
%    ind    :  index to selected elements

switch path
  case 'semi'
    data = struct( 'r', r, 'z', z, 'k0', k0, 'k', k );
    ind = 1 : numel( r );
  case 'real'
    ind = z >= obj.ratio * r;
    data = struct( 'r', r( ind ), 'z', z( ind ), 'k0', k0, 'k', k );
  case 'imag'
    ind = z < obj.ratio * r;
    data = struct( 'r', r( ind ), 'z', z( ind ), 'k0', k0, 'k', k );
end
