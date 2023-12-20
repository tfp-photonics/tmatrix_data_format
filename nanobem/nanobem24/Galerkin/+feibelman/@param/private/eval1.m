function [ d1, d2 ] = eval1( obj, k0 )
%  EVAL - Evaluate Feibelman parameters.
%
%  Usage for obj = feibelman.param :
%    obj = eval1( obj, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum

%  Feibelman function uses eV, internal functions use k0 (1/nm)
eV2nm = 1 / 8.0655477e-4;
[ d1, d2 ] = obj.fun( eV2nm / ( 2 * pi ) * k0 );

%  expand to full size, if needed
n = numel( obj.tau );
if size( d1, 1 ) ~= n
  d1 = repelem( d1, n, 1 );
  d2 = repelem( d2, n, 1 );
end
