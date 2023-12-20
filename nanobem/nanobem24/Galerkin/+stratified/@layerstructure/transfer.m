function m = transfer( obj, k0, kpar, varargin )
%  TRANSFER - Transfer matrix for given interface.
%
%  Usage for obj = stratified.layerstructure :
%    m = transfer( obj, k0, kpar, i1, i2 )
%  Input
%    k0     :  wavenumber of light in vacuum
%    kpar   :  parallel momenta
%    i1     :  first material index
%    i2     :  second material index
%  Output
%    m.te   :  transfer matrix for TE polarization 
%    m.tm   :  transfer matrix for TM polarization 

%  reflection and transmission coefficients
[ r, t ] = fresnel( obj, k0, kpar, varargin{ : } );
%  transfer matrix between layers, Hohenester Eq. (8.30)
te = cat( 3, 1 ./ t.te, r.te ./ t.te, r.te ./ t.te, 1 ./ t.te );
tm = cat( 3, 1 ./ t.tm, r.tm ./ t.tm, r.tm ./ t.tm, 1 ./ t.tm );

m.te = reshape( permute( te, [ 3, 1, 2 ]  ), 2, 2, [] );
m.tm = reshape( permute( tm, [ 3, 1, 2 ]  ), 2, 2, [] );
