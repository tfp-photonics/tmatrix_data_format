function plot( obj, varargin )
%  PLOT - Plot tabulated Green function.
    
%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'mode', 1 );
addParameter( p, 'iz', 1 );
%  parse input
parse( p, varargin{ : } );

if ~strcmp( p.UsingDefaults, 'mode' ), varargin = varargin( 2 : end );  end

switch p.Results.mode
  case 1
    z = obj.ztab1;
    [ rr, zz ] = ndgrid( obj.rtab, z );
    for name = fieldnames( obj.ytab ) .'
      ytab.( name{ 1 } ) = griddedInterpolant( rr, zz, squeeze(  ...
        obj.ytab.( name{ 1 } ).Values( :, :, p.Results.iz ) ), 'cubic', 'none' );
    end
  case 2
    z = obj.ztab2;
    [ rr, zz ] = ndgrid( obj.rtab, z );
    for name = fieldnames( obj.ytab ) .'
      ytab.( name{ 1 } ) = griddedInterpolant( rr, zz, squeeze(  ...
        obj.ytab.( name{ 1 } ).Values( :, p.Results.iz, : ) ), 'cubic', 'none' );
    end    
end

if isempty( obj.ytab )
  warning( 'stratified.tab.inter YTAB not initialized' );
else
  stratified.plot( obj.rtab, z, ytab, varargin{ : } );
end
