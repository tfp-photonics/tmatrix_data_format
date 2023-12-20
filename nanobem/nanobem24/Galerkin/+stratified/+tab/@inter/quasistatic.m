function yout = quasistatic( obj, pos1, pos2, k0, varargin )
%  QUASISTATIC - Quasistatic contribution for Green function elements.
%
%  Usage for obj = stratified.tab.inter :
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

if abs( obj.i1 - obj.i2 ) ~= 1
  yout = [];
else
  %  convert to coordinates for stratified system
  [ r, Z, R, siz ] = cart2strat( obj, pos1, pos2, varargin{ : } );
  r( r < 1e-10 ) = 1e-10;
  %  propagation direction of transmitted waves
  if obj.i1 < obj.i2
    dir = - 1;   
  else
    dir = + 1;
  end  
  
  %  first and second derivative of Green function along z
  y.z1 = dir / ( 4 * pi ) * log( Z + R );
  y.zz = - 1 ./ ( 4 * pi * R );
  %  surface Green function, Chew (6.28)
  y.s  = 1 ./ ( 4 * pi * R );
  %  components for double layer potential
  y.r   = - ( R - Z ) ./ ( 4 * pi * r );
  y.rz1 = dir * ( 1 - Z ./ R ) ./ ( 4 * pi * r ); 
  %  derivatives wrt Z2
  [ y.z2, y.rz2 ] = deal( - y.z1, - y.rz1 ); 

  switch p.Results.refl
    case 0
      %  reshape Green functions
      for name = [ "z1", "z2", "zz", "s", "r", "rz1", "rz2" ]
        yout.( name ) = reshape( y.( name{ 1 } ), siz );
      end
    case 1
      %  transmission coefficient for quasistatic approximation
      tq = fresnelquasistatic( obj, k0 );
      %  reshape Green functions and multiply with prefactor
      for name1 = [ "te", "tm" ]
      for name2 = [ "z1", "z2", "zz", "s", "r", "rz1", "rz2" ]
        if abs( tq.( name1 ) ) > 1e-10
          yout.( name1 + name2 ) = reshape( tq.( name1 ) * y.( name2 ), siz );
        end
      end
      end
  end 
  %  deal with empty structure
  if ~exist( 'yout', 'var' ),  yout = [];  end
end
