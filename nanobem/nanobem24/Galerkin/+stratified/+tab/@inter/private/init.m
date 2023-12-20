function obj = init( obj, i1, i2, rtab, ztab1, ztab2, varargin )
%  INIT - Initialize tabulated Green functions connecting different layers.
%
%  Usage for obj = stratified.tab.inter :
%    obj = init( obj, i1, i2, rtab, ztab1, ztab2, PropertyPairs )
%  Input
%    i1       :  medium index for observation points
%    i2       :  medium index for source points
%    rtab     :  radii for tabulated Green functions
%    ztab1    :  z-values observation points of tabulated Green functions 
%    ztab2    :  z-values source points of tabulated Green functions 

[ obj.i1, obj.i2, obj.rtab, obj.ztab1, obj.ztab2 ] =  ...
               deal( i1, i2, rtab, unique( ztab1 ), unique( ztab2 ) );
             
%  deal with single z values, griddedInterpolant with cubic interpolation
%  requires at least three values
eta = 0.1;
if isscalar( obj.ztab1 ),  obj.ztab1 = obj.ztab1 + [ - eta, 0, eta ];  end
if isscalar( obj.ztab2 ),  obj.ztab2 = obj.ztab2 + [ - eta, 0, eta ];  end
