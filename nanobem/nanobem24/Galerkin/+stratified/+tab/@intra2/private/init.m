function obj = init( obj, i1, rtab, ztab1, ztab2, varargin )
%  INIT - Initialize tabulated Green functions for inside layer.
%
%  Usage for obj = stratified.tab.intra2 :
%    obj = init( obj, layer, i1, rtab, ztab1, ztab2 )
%  Input
%    layer    :  layer structure
%    i1       :  layer index
%    rtab     :  radii for tabulated Green functions
%    ztab1    :  z-values for tabulated Green function, z1-z2
%    ztab2    :  z-values for tabulated Green function, z1+z2
      
[ obj.i1, obj.rtab, obj.ztab1, obj.ztab2 ] = deal( i1, rtab, ztab1, ztab2 );
