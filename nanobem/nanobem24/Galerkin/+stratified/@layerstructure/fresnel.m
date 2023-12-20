function [ r, t ] = fresnel( obj, k0, kpar, i1, varargin )
%  FRESNEL - Fresnel coefficients for given interface.
%
%  Usage for obj = stratified.layerstructure :
%    [ r, t ] = fresnel( obj, k0, kpar, i1, i2, PropertyPairs )
%  Input
%    k0     :  wavenumber of light in vacuum
%    kpar   :  parallel momenta
%    i1     :  first  material index
%    i2     :  second material index
%  PropertyName
%    efield :  TM transmission coefficient for electric or magnetic field
%  Output
%    r      :  reflection coefficients
%    t      :  transmission coefficients

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'i2', i1 + 1 );
addParameter( p, 'efield', 0 );
%  parse input
parse( p, varargin{ : } );

%  material properties
obj = eval( obj, k0 );
i2 = p.Results.i2;
[ eps1, mu1 ] = deal( obj.data.eps( i1 ), obj.data.mu( i1 ) );
[ eps2, mu2 ] = deal( obj.data.eps( i2 ), obj.data.mu( i2 ) );
%  wavenumber
kz1 = stratified.zsqrt( eps1 * mu1 * k0 ^ 2 - kpar .^ 2 );
kz2 = stratified.zsqrt( eps2 * mu2 * k0 ^ 2 - kpar .^ 2 );
%  impediancies
Z1 = zsqrt( mu1 / eps1 );
Z2 = zsqrt( mu2 / eps2 );

%  reflection and transmission coefficients for TE polarization
r.te = ( mu2 * kz1 - mu1 * kz2 ) ./ ( mu2 * kz1 + mu1 * kz2 );
t.te =           2 * mu2 * kz1   ./ ( mu2 * kz1 + mu1 * kz2 );
%  reflection and transmission coefficients for TM polarization
r.tm = ( eps2 * kz1 - eps1 * kz2 ) ./ ( eps2 * kz1 + eps1 * kz2 );
t.tm =            2 * eps2 * kz1   ./ ( eps2 * kz1 + eps1 * kz2 );
%  TM coefficient for electric or magnetic field
if p.Results.efield,  t.tm = t.tm * Z2 / Z1;  end

