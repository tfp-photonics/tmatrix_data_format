function yout = quasistatic( obj, pos1, pos2, k0, varargin )
%  QUASISTATIC - Quasistatic contribution for Green function elements.
%
%  Usage for obj = stratified.tab.intra2 :
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
[ r, Z1, Z2, R1, R2, siz ] = cart2strat( obj, pos1, pos2, varargin{ : } );
r( r < 1e-10 ) = 1e-10;

%  first and second derivative of Green function along z
y.z1  =   1 / ( 4 * pi ) * log( Z1 + R1 );
y.z2  = - 1 / ( 4 * pi ) * log( Z2 + R2 );
y.zz1 = 1 ./ ( 4 * pi * R1 );
y.zz2 = 1 ./ ( 4 * pi * R2 );
%  surface Green function, Chew (6.28)
y.s1  = 1 ./ ( 4 * pi * R1 );
y.s2  = 1 ./ ( 4 * pi * R2 );
%  components for double layer potential
y.r1  = - ( R1 - Z1 ) ./ ( 4 * pi * r );
y.r2  = - ( R2 - Z2 ) ./ ( 4 * pi * r );
y.rz1 =   ( 1 - Z1 ./ R1 ) ./ ( 4 * pi * r );
y.rz2 = - ( 1 - Z2 ./ R2 ) ./ ( 4 * pi * r );

switch p.Results.refl
  case 0
    %  reshape Green functions
    for name = [ "z", "zz", "s", "r", "rz" ]
      yout( 1 ).( name ) = reshape( y.( name + "1" ), siz );
      yout( 2 ).( name ) = reshape( y.( name + "2" ), siz );
    end
  case 1
    %  reflection coefficient for quasistatic approximation
    rq = fresnelquasistatic( obj, k0 );
    %  reshape Green functions and multiply with prefactor
    for name1 = [ "te", "tm" ]
    for name2 = [ "z", "zz", "s", "r", "rz" ]
      if abs( rq( 1 ).( name1 ) ) + abs( rq( 2 ).( name1 ) ) > 1e-10
        yout.( name1 + name2 ) = reshape(  ...
          rq( 1 ).( name1 ) * y.( name2 + "1" ) +  ...
          rq( 2 ).( name1 ) * y.( name2 + "2" ), siz );
      end
    end
    end
    %  deal with empty structure
    if ~exist( 'yout', 'var' ),  yout = [];  end    
end
