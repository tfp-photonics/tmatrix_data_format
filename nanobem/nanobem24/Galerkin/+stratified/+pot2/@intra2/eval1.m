function data = eval1( obj, green, k0, varargin )
%  EVAL1 - Evaluate w/o refinement for all boundary element pairs.
%
%  Usage for obj = stratified.pot2.intra2 :
%    data = eval2( obj, green, k0, PropertyPairs )
%  Input
%    green  :  tabulated Green function object
%    k0     :  wavenumber of light in vacuum
%  Output
%    data   :  single and double layer potential

%  integration points, weights and shape functions
[ i1, i2, q2, a2, k ] = deal( 1, 2, 3, 4, 5 );
data = oint1( obj, k0, 'ind', [ i1, i2, q2, a2, k ], varargin{ : } );
%  evaluate Green function elements 
[ y1, y2 ] = interp1( obj, green,  ...
    data.pos1, data.pos2, k0, 'ind', [ i1, i2, q2, k ], 'smooth', 1 );  
[ z1, z2, data1 ] = refine2( obj, data, k0, 'ind', [ i1, i2, q2, k ] );
%  add smooth and quasistatic elements
if ~isempty( z1 )
  for name = convertCharsToStrings( fieldnames( z1 ) ) .'
    y1.( name ) = y1.( name ) + z1.( name );
    y2.( name ) = y2.( name ) + z2.( name );
  end
end  

%  unit vectors in radial, azimuthal and z-direction
[ er, et, ez ] = deal( data.er, data.et, data.ez );
%  shape functions
[ f, fp ] = deal( data.f, data.fp );
fr = dot( er, f, k );
ft = dot( et, f, k );
fz = dot( ez, f, k );

%  auxiliary quantity
x = ( y1.tez  + y1.tmz ) * ez * fp - ( y2.tmz + y2.tez ) * fz;
%  electric single layer potential
SL1 = - x +                                        ...
  ( - y2.tmzz / data.k1 ^ 2 + y2.te ) * fp +       ...
  (   y1.tm * data.k1 ^ 2 - y1.tezz ) * ez * fz +  ... 
      y1.tes * ( fr * er + et * ft );
%  magnetic single layer potential 
SL2 = - x +                                        ...
  ( - y2.tezz / data.k1 ^ 2 + y2.tm ) * fp +       ...
  (   y1.te * data.k1 ^ 2 - y1.tmzz ) * ez * fz +  ... 
      y1.tms * ( er * fr + et * ft );
      
%  electric double layer potential
DL1 = y1.tmrz * et * fr ./ data.r -                         ...
      y2.tmrz *      ft - y1.tmr * ez * ft * data.k1 ^ 2 -  ...
      y1.terz * et * fp + y1.ter * et * fz * data.k1 ^ 2;
%  magnetic double layer potential
DL2 = y1.terz * et * fr ./ data.r -                         ...
      y2.terz *      ft - y1.ter * ez * ft * data.k1 ^ 2 -  ...
      y1.tmrz * et * fp + y1.tmr * et * fz * data.k1 ^ 2;    

%  material parameters
[ mu1, eps1 ] = deal( data.mu1, data.eps1 );
%  perform integration
data = struct;
data.SL1 = double( sum( SL1, q2 ), [ i1, i2, a2, k ] );
data.SL2 = double( sum( SL2, q2 ), [ i1, i2, a2, k ] );
data.DL1 = double( sum( DL1, q2 ), [ i1, i2, a2, k ] );
data.DL2 = double( sum( DL2, q2 ), [ i1, i2, a2, k ] );
%  add refined quasistatic elements
if ~isempty( data1 )
  data.SL1 = data.SL1 + data1.SL1;  data.SL2 = data.SL2 + data1.SL2;
  data.DL1 = data.DL1 + data1.DL1;  data.DL2 = data.DL2 + data1.DL2;
end
%  multiply with prefactors
data.SL1 = data.SL1 * mu1;
data.SL2 = data.SL2 * eps1;
