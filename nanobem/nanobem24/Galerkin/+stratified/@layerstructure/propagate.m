function p = propagate( obj, k0, kpar, i1 )
%  PROPAGATE - Propagation matrix through given layer.
%
%  Usage for obj = stratified.layerstructure :
%    p = propagate( obj, k0, kpar, i1 )
%  Input
%    k0     :  wavenumber of light in vacuum
%    kpar   :  parallel momenta
%    i1     :  index to layer
%  Output
%    p      :  propagation matrix

%  wavenumber and thickness of layer
kz = subsref( obj, substruct( '.', 'kz', '()', { k0, kpar, i1 } ) );
d  = obj.d( i1 - 1 );

%  propagation matrix
%    For thick layer or wavenumbers with large imaginary parts one can
%    encounter problems with too large arguments.  Rather than
%    introducing a description in terms of the scattering matrix we
%    here proceed in a more simple fashion.
p = cat( 3, exp( - 1i * kz * d ), 0 * kz, 0 * kz, exp( 1i * kz * d ) );
p = reshape( permute( p, [ 3, 1, 2 ] ), 2, 2, [] );
p( isinf( p ) ) = 1e30;
