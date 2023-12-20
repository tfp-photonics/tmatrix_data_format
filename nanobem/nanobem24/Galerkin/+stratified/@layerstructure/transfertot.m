function mtot = transfertot( obj, k0, kpar, varargin )
%  TRANSFERTOT - Total transfer matrix for given layer structure.
%
%  Usage for obj = stratified.layerstructure :
%    mtot = transfertot( obj, k0, kpar )
%  Input
%    k0     :  wavenumber of light in vacuum
%    kpar   :  parallel momenta
%  Output
%    mtot   :  total transfer matrix

%  transfer matrix at first interface, Hohenester Eq. (8.32)
mtot = transfer( obj, k0, kpar, 1, varargin{ : } );
%  loop through layer structure
for it = 2 : obj.n
  %  propagation and transfer matrices
  p = propagate( obj, k0, kpar, it ); 
  m = transfer(  obj, k0, kpar, it, varargin{ : } );
  %  multiply transfer matrices
  mtot.te = fun( mtot.te, fun( p, m.te ) );
  mtot.tm = fun( mtot.tm, fun( p, m.tm ) );
end


function z = fun( x, y )
%  Multi-dimensional multiplication function.
switch ndims( x )
  case 2
    z = x * y;
  otherwise
    z = cat( 1,  ...
      x( 1, 1, : ) .* y( 1, 1, : ) + x( 1, 2, : ) .* y( 2, 1, : ),  ...
      x( 2, 1, : ) .* y( 1, 1, : ) + x( 2, 2, : ) .* y( 2, 1, : ),  ...
      x( 1, 1, : ) .* y( 1, 2, : ) + x( 1, 2, : ) .* y( 2, 2, : ),  ...
      x( 2, 1, : ) .* y( 1, 2, : ) + x( 2, 2, : ) .* y( 2, 2, : ) );
    z = reshape( z, 2, 2, [] );
end