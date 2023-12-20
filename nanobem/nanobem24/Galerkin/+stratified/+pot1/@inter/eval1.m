function data = eval1( obj, green, k0, varargin )
%  EVAL1 - Evaluate w/o refinement for all boundary element pairs.
%
%  Usage for obj = stratified.pot1.intra1 :
%    data = eval1( obj, green, k0, PropertyPairs )
%  Input
%    green  :  tabulated Green function object
%    k0     :  wavenumber of light in vacuum
%  Output
%    data   :  single and double layer potential

%  boundary integrators and quadrature points
[ oint, data ] = oint1( obj, k0 );
%  evaluate Green function elements 
y = interp( green, data.pos1, data.pos2, k0, 'smooth', 1 );
[ y1, data1 ] = refine2( obj, data, k0 );
%  add smooth and quasistatic elements
if ~isempty( y1 )
  for name = intersect( convertCharsToStrings( fieldnames( y  ) ),  ...
                        convertCharsToStrings( fieldnames( y1 ) ) ) .'
    y.( name ) = y.( name ) + y1.( name );
  end
end
%  material parameters
k1 = data.k1;  eps1 = data.eps1;  mu1 = data.mu1;
k2 = data.k2;  eps2 = data.eps2;  mu2 = data.mu2;
%  ratio of impedances
z = sqrt( mu1 / eps1 ) / sqrt( mu2 / eps2 );
data = struct;

for a1 = 1 : obj.pt1.npoly
for a2 = 1 : obj.pt2.npoly
  %  electric single layer potential, Chew (6.28)
  data.SL1( :, a1, :, a2 ) =  ...
    oint.pp( a1, a2, y.tmzz * z / k1 / k2 - y.te   ) -  ...
    oint.zp( a1, a2, y.tmz2 * z * k1 / k2 + y.tez1 ) -  ...
    oint.pz( a1, a2, y.tmz1 * z / k1 * k2 + y.tez2 ) +  ...
    oint.zz( a1, a2, y.tm   * z * k1 * k2 - y.tezz ) +  ...
    oint.ss( a1, a2, y.tes );
  %  magnetic single layer potential
  data.SL2( :, a1, :, a2 ) =  ...
    oint.pp( a1, a2, y.tezz / z / k1 / k2 - y.tm   ) -  ...
    oint.zp( a1, a2, y.tez2 / z * k1 / k2 + y.tmz1 ) -  ...
    oint.pz( a1, a2, y.tez1 / z / k1 * k2 + y.tmz2 ) +  ...
    oint.zz( a1, a2, y.te   / z * k1 * k2 - y.tmzz ) +  ...
    oint.ss( a1, a2, y.tms );

  %  electric double layer potential
  data.DL1( :, a1, :, a2 ) =  ...  
    oint.pr( a1, a2, y.tmrz1 * z ) / k1 * k2 - oint.zr( a1, a2, y.tmr * z ) * k1 * k2 -  ...
    oint.rp( a1, a2, y.terz2     )           + oint.rz( a1, a2, y.ter     ) * k2 * k2;   
  %  magnetic double layer potential, Chew (6.41), (6.43)
  data.DL2( :, a1, :, a2 ) =  ...  
    oint.pr( a1, a2, y.terz1 / z ) / k1 * k2 - oint.zr( a1, a2, y.ter / z ) * k1 * k2 -  ...
    oint.rp( a1, a2, y.tmrz2     )           + oint.rz( a1, a2, y.tmr     ) * k2 * k2; 
end
end

%  add refined quasistatic elements
if ~isempty( data1 )
  data.SL1 = data.SL1 + data1.SL1;  data.SL2 = data.SL2 + data1.SL2;
  data.DL1 = data.DL1 + data1.DL1;  data.DL2 = data.DL2 + data1.DL2;
end
%  multiply with prefactors
data.SL1 = data.SL1 * mu2;
data.SL2 = data.SL2 * eps2; 
