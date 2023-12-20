function [ d, ind ] = mindist( obj, pos )
%  MINDIST - Minimal z-distance to interfaces.
%
%  Usage for obj = stratified.layerstructure :
%    [ d, ind ] = mindist( obj, pos )
%  Input
%    pos    :  positions
%  Output
%    d      :  minimal z-distance to interface
%    ind    :  interface index

%  distances from positions to layer interfaces
d = abs( bsxfun( @minus, pos( :, 3 ), reshape( obj.z, 1, [] ) ) );
%  find mimimum
[ d, ind ] = min( d, [],  2 );
