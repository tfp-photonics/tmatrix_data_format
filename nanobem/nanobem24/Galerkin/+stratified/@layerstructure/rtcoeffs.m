function [ r, t ] = rtcoeffs( obj, k0, kpar, varargin )
%  RTCOEFFS - Reflection and transmission coefficients for layer structure.
%
%  Usage for obj = stratified.layerstructure :
%    [ r, t ] = rtcoeffs( obj, k0, kpar, PropertyName )
%  Input
%    k0     :  wavenumber of light in vacuum
%    kpar   :  parallel momenta
%  PropertyName
%    dir    :  'up' for upgoing and 'down' for downgoing wave
%  Output
%    r      :  reflection coefficients for layer structure
%    t      :  transmission coefficients for layer structure

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'dir', 'up' );
%  parse input
parse( p, varargin{ : } );

%  transfer matrix for layer structure
mtot = transfertot( obj, k0, kpar );
%  reflection and transmission coefficients
switch p.Results.dir
  case 'up'
    %  upgoing wave, Hohenester Eq. (8.34)
    r.te = mtot.te( 2, 1, : ) ./ mtot.te( 1, 1, : );
    t.te =                  1 ./ mtot.te( 1, 1, : );

    r.tm = mtot.tm( 2, 1, : ) ./ mtot.tm( 1, 1, : );
    t.tm =                  1 ./ mtot.tm( 1, 1, : );
    
  case { 'dn', 'down' }
    %  downgoing wave, Hohenester Eq. (8.36)
    r.te = - mtot.te( 1, 2, : ) ./ mtot.te( 1, 1, : );
    t.te =   mtot.te( 2, 2, : )  + mtot.te( 2, 1, : ) .* r.te;
    
    r.tm = - mtot.tm( 1, 2, : ) ./ mtot.tm( 1, 1, : );
    t.tm =   mtot.tm( 2, 2, : )  + mtot.tm( 2, 1, : ) .* r.tm; 
end

%  reshape output
siz = size( kpar );
r = struct( 'te', reshape( r.te, siz ), 'tm', reshape( r.tm, siz ) );
t = struct( 'te', reshape( t.te, siz ), 'tm', reshape( t.tm, siz ) );
