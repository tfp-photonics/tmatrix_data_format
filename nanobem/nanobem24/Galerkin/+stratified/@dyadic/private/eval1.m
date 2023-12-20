function data = eval1( tab, pos1, pos2, k0, varargin )
%  EVAL1 - Evaluate dyadic intralayer Green function.
%
%  Usage for obj = stratified.dyadic : 
%    data = eval1( obj, pos1, pos2, k0, PropertyPairs );
%  Input
%    tab      :  tabulated Green function
%    pos1     :  observation points in single layer
%    pos2     :  source points in the same layer
%    k0       :  wavenumber of light in vacuum
%  Output
%    data     :  single and double layer potentials

%  Green functions and auxiliary data as tensor object
[ i1, i2, k1, k2 ] = deal( 1, 2, 3, 4 );
[ y1, y2, data ] = tensorize( tab, pos1, pos2, k0,  ...
                             varargin{ : }, 'ind', [ i1, i2, k1, k2 ] );
                              
%  auxiliary quantity
x = - ( y2( 2 ).tez + y2( 2 ).tmz ) * data.ez1 -  ...
      ( y2( 1 ).tmz + y2( 1 ).tez ) * data.ez2;
%  electric single layer potential
SL1 = - x +  ...
  ( y2( 3 ).tmzz / data.k1 ^ 2 - y2( 3 ).te ) +  ...
  ( y1.tm * data.k1 ^ 2 - y1.tezz ) * data.ez1 * data.ez2 +  ...
    y1.tes * (  data.ex1 * data.ex2 + data.ey1 * data.ey2 );
%  magnetic single layer potential
SL2 = - x +  ...
  ( y2( 3 ).tezz / data.k1 ^ 2 - y2( 3 ).tm ) +  ...
  ( y1.te * data.k1 ^ 2 - y1.tmzz ) * data.ez1 * data.ez2 +  ...
    y1.tms * (  data.ex1 * data.ex2 + data.ey1 * data.ey2 );  

%  electric double layer potential 
DL1 = cross( data.ez2, y2( 3 ).tmz, k2 ) - y1.tmr * data.ez1 * data.et2 * data.k1 ^ 2 +  ...
      cross( data.ez1, y2( 3 ).tez, k1 ) + y1.ter * data.et1 * data.ez2 * data.k1 ^ 2;
%  magnetic double layer potential    
DL2 = cross( data.ez2, y2( 3 ).tez, k2 ) - y1.ter * data.ez1 * data.et2 * data.k1 ^ 2 +  ...
      cross( data.ez1, y2( 3 ).tmz, k1 ) + y1.tmr * data.et1 * data.ez2 * data.k1 ^ 2;
      
%  material parameters
[ mu1, eps1 ] = deal( data.mu1, data.eps1 );
%  save output, multiply with prefactors
data = struct;
data.SL1 = double( SL1, [ i1, k1, i2, k2 ] ) * mu1;
data.SL2 = double( SL2, [ i1, k1, i2, k2 ] ) * eps1;
data.DL1 = double( DL1, [ i1, k1, i2, k2 ] );
data.DL2 = double( DL2, [ i1, k1, i2, k2 ] );
