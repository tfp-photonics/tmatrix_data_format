function data = eval1( obj, green, k0, varargin )
%  EVAL1 - Evaluate w/o refinement for all boundary element pairs.
%
%  Usage for obj = stratified.pot1.inter :
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
%  material parameters
k1 = data.k1;  eps1 = data.eps1;  mu1 = data.mu1;
k2 = data.k2;  eps2 = data.eps2;  mu2 = data.mu2; 
%  ratio of impedances
z = sqrt( mu1 / eps1 ) / sqrt( mu2 / eps2 );

%  electric single layer potential
SL1 = ( - y2.tmzz * z / k1 / k2 + y2.te   ) *      fp -  ...
      (   y1.tmz2 * z * k1 / k2 + y1.tez1 ) * ez * fp -  ...
      ( - y2.tmz1 * z / k1 * k2 - y2.tez2 ) *      fz +  ...
      (   y1.tm   * z * k1 * k2 - y1.tezz ) * ez * fz +  ...
          y1.tes  * ( er * fr + et * ft );
%  magnetic single layer potential 
SL2 = ( - y2.tezz / z / k1 / k2 + y2.tm   ) *      fp -  ...
      (   y1.tez2 / z * k1 / k2 + y1.tmz1 ) * ez * fp -  ...
      ( - y2.tez1 / z / k1 * k2 - y2.tmz2 )      * fz +  ...
      (   y1.te   / z * k1 * k2 - y1.tmzz ) * ez * fz +  ...
          y1.tms  * ( er * fr + et * ft );

%  electric double layer potential
DL1 = y1.tmrz1 * z * et * fr * k2 / k1 ./ data.r -                         ...
      y2.tmrz1 * z      * ft * k2 / k1 - y1.tmr * z * ez * ft * k1 * k2 -  ...
      y1.terz2     * et * fp           + y1.ter     * et * fz * k2 * k2;             
%  magnetic double layer potential
DL2 = y1.terz1 / z * et * fr * k2 / k1 ./ data.r -                         ...
      y2.terz1 / z      * ft * k2 / k1 - y1.ter / z * ez * ft * k1 * k2 -  ...
      y1.tmrz2     * et * fp           + y1.tmr     * et * fz * k2 * k2;

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
data.SL1 = data.SL1 * mu2;
data.SL2 = data.SL2 * eps2; 
