function [ y1, y2, data ] = refine2( obj, data1, k0, varargin )
%  REFINE2 - Handle refined quasistatic terms.
%
%  Usage for obj = stratified.pot2.inter :
%    [ y1, y2, data ] = refine2( obj, data1, k0, PropertyPairs )
%  Input
%    data1  :  positions and material properties
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    ind    :  tensor indices [i1,i2,q2,k]
%  Output
%    y1     :  quasistatic Green function elements
%    y2     :  derivatives of quasistatic Green function elements
%    data   :  refined single and double layer potential terms

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ind', 1 : 4 );
%  parse input
parse( p, varargin{ : } );

%  evaluate quasistatic Green function elements and derivatives
green = stratified.tab.inter( obj.layer, obj.i1, obj.i2, [], [], [] );
[ y1, y2 ] = interp1( obj, green,  ...
   data1.pos1, data1.pos2, k0, 'stat', 1, 'refl', 1, varargin{ : } );  

%  deal with refined Green function elements
if isempty( obj.yout ) || isempty( y1 )
  data = [];
else
  %  set refined Green function elements to zero
  yout = obj.yout;
  [ y1, siz ] = fzero( obj, y1, p.Results.ind( 1 : 3 ), yout.i1, yout.i2 );
  [ y2, ~   ] = fzero( obj, y2, p.Results.ind( 1 : 4 ), yout.i1, yout.i2 );
  %  wavenumbers 
  k1 = data1.k1;
  k2 = data1.k2;
  %  ratio of impedances
  z = sqrt( data1.mu1 / data1.eps1 ) / sqrt( data1.mu2 / data1.eps2 );
  
  %  quasistatic transmission coefficients
  tq = fresnelquasistatic( green, k0 );
  %  refined electric single layer elements
  SL1 = - tq.tm * z * yout.pp  / k1 / k2 -                     ...
          tq.tm * z * yout.zp2 * k1 / k2 - tq.te * yout.zp1 +  ...
          tq.tm * z * yout.pz1 / k1 * k2 + tq.te * yout.pz2 -  ...
          tq.te * yout.zz + tq.te * yout.ss;
  %  refined electric single layer elements 
  SL2 = - tq.te / z * yout.pp  / k1 / k2 -                     ...
          tq.te / z * yout.zp2 * k1 / k2 - tq.tm * yout.zp1 +  ...
          tq.te / z * yout.pz1 / k1 * k2 + tq.tm * yout.pz2 -  ...
          tq.tm * yout.zz + tq.tm * yout.ss;
           
  %  refined electric double layer potentials
  DL1 = tq.tm * z * ( yout.pr * k2 / k1 - yout.zr * k1 * k2 ) -  ...
        tq.te     * ( yout.rp           - yout.rz * k2 * k2 );
  %  refined magnetic double layer potentials
  DL2 = tq.te / z * ( yout.pr * k2 / k1 - yout.zr * k1 * k2 ) -  ...
        tq.tm     * ( yout.rp           - yout.rz * k2 * k2 );  
  
  %  allocate output
  [ data.SL1, data.SL2, data.DL1, data.DL2 ] =  ...
                   deal( zeros( siz( 1 ) * siz( 2 ), size( SL1, 2 ), 3 ) );
  %  set refined elements
  ind = sub2ind( siz( [ 1, 2 ] ), yout.i1, yout.i2 );
  data.SL1( ind, :, : ) = SL1;  
  data.SL2( ind, :, : ) = SL2;  
  data.DL1( ind, :, : ) = DL1;  
  data.DL2( ind, :, : ) = DL2;  
  %  reshape output
  for name = [ "SL1", "SL2", "DL1", "DL2" ]
    data.( name ) = reshape( data.( name ), siz( 1 ), siz( 2 ), [], 3 );
  end
end
