function plot( obj, varargin )
%  PLOT - Plot tabulated Green function.
  
%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'z', 1 );
%  parse input
parse( p, varargin{ : } );

if isempty( p.UsingDefaults ), varargin = varargin( 2 : end );  end

if isempty( obj.ytab1 )
  warning( 'stratified.tab.intra2 YTAB not initialized' );
else
  switch p.Results.z
    case 1
      stratified.plot( obj.rtab, obj.ztab1, obj.ytab1, varargin{ : } );
    otherwise
      stratified.plot( obj.rtab, obj.ztab2, obj.ytab2, varargin{ : } );
  end
end
