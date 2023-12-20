function yout = quasistatic( obj, pos1, pos2, k0, varargin )
%  QUASISTATIC - Quasistatic contribution for Green function elements.
%
%  Usage for obj = stratified.tab.intra1 :
%    yout = quasistatic( obj, pos1, pos2, k0, varargin )
%  Input
%    pos1   :  observation points
%    pos2   :  source points
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    rows   :  position arrays of same size
%    refl   :  multiply with quasistatic reflection coefficients
%  Output
%    yout  :  quasistatic Green function elements

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'refl', 0 );
%  parse input
parse( p, varargin{ : } );

%  convert to coordinates for stratified system
[ r, Z, R, siz ] = cart2strat( obj, pos1, pos2, varargin{ : } );
r( r < 1e-10 ) = 1e-10;
%  propagation direction of reflected waves
if obj.i1 == 1
  dir = - 1;  
else
  dir = + 1;
end

%  first and second derivative of Green function along z
y.z  = dir / ( 4 * pi ) * log( Z + R );
y.zz = 1 ./ ( 4 * pi * R );
%  surface Green function, Chew (6.28)
y.s  = 1 ./ ( 4 * pi * R );
%  components for double layer potential
y.r  = - ( R - Z ) ./ ( 4 * pi * r );
y.rz = dir * ( 1 - Z ./ R ) ./ ( 4 * pi * r );

switch p.Results.refl
  case 0
    %  reshape Green functions
    for name = [ "z", "zz", "s", "r", "rz" ]
      yout.( name ) = reshape( y.( name{ 1 } ), siz );
    end
  case 1
    %  reflection coefficient for quasistatic approximation
    rq = fresnelquasistatic( obj, k0 );
    %  reshape Green functions and multiply with prefactor
    for name1 = [ "te", "tm" ]
    for name2 = [ "z", "zz", "s", "r", "rz" ]
      if abs( rq.( name1 ) ) > 1e-10
        yout.( name1 + name2 ) = reshape( rq.( name1 ) * y.( name2 ), siz );
      end
    end
    end
    %  deal with empty structure
    if ~exist( 'yout', 'var' ),  yout = [];  end
end
