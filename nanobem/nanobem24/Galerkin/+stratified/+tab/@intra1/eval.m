function y = eval( obj, pos1, pos2, k0, varargin )
%  EVAL - Evaluate Green function elements.
%
%  Usage for obj = stratified.tab.intra1 :
%    y = eval( obj, pos1, pos2, k0, varargin )
%  Input
%    pos1     :  observation points
%    pos2     :  source points
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    name     :  names of reflected Green functions
%    rows     :  positions of same size or not
%    smooth   :  w/o or with quasistatic contribution

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'smooth', 0 );
%  parse input
parse( p, varargin{ : } );

%  fill Green function tables ?
if isempty( obj.k0 ) || k0 ~= obj.k0,  obj = fill( obj, k0 );  end

%  convert to coordinates for stratified system
[ r, Z, ~, siz ] = cart2strat( obj, pos1, pos2, varargin{ : } );
r( r < min( obj.rtab ) ) = min( obj.rtab );
%  perform interpolation
for name = convertCharsToStrings( fieldnames( obj.ytab ) ) .'
  y.( name ) = reshape( obj.ytab.( name )( r, Z ), siz );
end

%  add quasistatic contribution ?
if ~p.Results.smooth
  %  compute quasistatic Green function
  y2 = quasistatic( obj, pos1, pos2, k0, varargin{ : }, 'refl', 1 );
  %  add to tabulated values
  if ~isempty( y2 )
    for name = convertCharsToStrings( fieldnames( y2 ) ) .'
      y.( name ) = y.( name ) + y2.( name ); 
    end
  end
end
