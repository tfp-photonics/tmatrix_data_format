function obj = init( obj, i1, rtab, ztab, varargin )
%  INIT - Initialize tabulated Green functions for upper or lower layer.
%
%  Usage for obj = stratified.tab.intra1 :
%    obj = init( obj, layer, i1, rtab, ztab )
%  Input
%    layer    :  layer structure
%    i1       :  layer index
%    rtab     :  radii for tabulated Green functions
%    ztab     :  z-values for tabulated Green functions 

[ obj.i1, obj.rtab, obj.ztab ] = deal( i1, rtab, ztab );
