function plot( obj, varargin )
%  PLOT - Plot tabulated Green function.
    
if isempty( obj.ytab )
  warning( 'stratified.tab.intra1 YTAB not initialized' );
else
  stratified.plot( obj.rtab, obj.ztab, obj.ytab, varargin{ : } );
end
