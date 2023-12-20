function [ y1, y2, data ] = refine2( obj, data1, k0, varargin )
%  REFINE2 - Handle refined quasistatic terms.
%
%  Usage for obj = stratified.pot2.intra1 :
%    [ y1, y2, data ] = refine2( obj, green, data1, k0, PropertyPairs )
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
green = stratified.tab.intra1( obj.layer, obj.i1, [], [] );
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
  
  %  quasistatic reflection coefficients
  rq = fresnelquasistatic( green, k0 );
  %  refined electric single layer elements
  SL1 = - ( rq.tm + rq.te ) * yout.x  -  ...
       rq.tm / data1.k1 ^ 2 * yout.pp -  ...
       rq.te * yout.zz + rq.te * yout.ss;
  %  refined magnetic single layer elements
  SL2 = - ( rq.tm + rq.te ) * yout.x  -  ...
       rq.te / data1.k1 ^ 2 * yout.pp -  ...
       rq.tm * yout.zz + rq.te * yout.ss;
           
  %  refined electric double layer potentials
  DL1 = rq.tm * yout.pr - rq.tm * yout.zr * data1.k1 ^ 2 -  ...
        rq.te * yout.rp + rq.te * yout.rz * data1.k1 ^ 2;
  %  refined magnetic double layer potentials
  DL2 = rq.te * yout.pr - rq.te * yout.zr * data1.k1 ^ 2 -  ...
        rq.tm * yout.rp + rq.tm * yout.rz * data1.k1 ^ 2;      
  
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
