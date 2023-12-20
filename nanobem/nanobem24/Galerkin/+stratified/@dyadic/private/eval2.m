function data = eval2( tab, pos1, pos2, k0, varargin )
%  EVAL2 - Evaluate dyadic interlayer Green function.
%
%  Usage for obj = stratified.dyadic : 
%    data = eval2( obj, pos1, pos2, k0, PropertyPairs );
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
                           
%  material parameters
eps1 = data.eps1;  mu1 = data.mu1;
eps2 = data.eps2;  mu2 = data.mu2;
%  ratio of impedances           
z = sqrt( mu1 / eps1 ) / sqrt( mu2 / eps2 );

%  electric single layer potential
SL1 = ( y2( 3 ).tmzz * z / data.k1 / data.k2 - y2( 3 ).te   )            +  ...
      ( y2( 2 ).tmz2 * z * data.k1 / data.k2 + y2( 2 ).tez1 ) * data.ez1 +  ...
      ( y2( 1 ).tmz1 * z / data.k1 * data.k2 + y2( 1 ).tez2 ) * data.ez2 +  ...
      ( y1.tm  * z * data.k1 * data.k2 - y1.tezz ) * data.ez1 * data.ez2 +  ...
        y1.tes * ( data.er1 * data.er2 + data.et1 * data.et2 );
%  magnetic single layer potential
SL2 = ( y2( 3 ).tezz / z / data.k1 / data.k2 - y2( 3 ).tm   )            +  ...
      ( y2( 2 ).tez2 / z * data.k1 / data.k2 + y2( 2 ).tmz1 ) * data.ez1 +  ...
      ( y2( 1 ).tez1 / z / data.k1 * data.k2 + y2( 1 ).tmz2 ) * data.ez2 +  ...
      ( y1.te  / z * data.k1 * data.k2 - y1.tmzz ) * data.ez1 * data.ez2 +  ...
        y1.tms * ( data.er1 * data.er2 + data.et1 * data.et2 );   
      
%  electric double layer potential
DL1 = cross( data.ez2, y2( 3 ).tmz1, k2 ) * z / data.k1 * data.k2 +  ...
      cross( data.ez1, y2( 3 ).tez2, k1 )                         -  ...
      y1.tmr * data.ez1 * data.et2        * z * data.k1 * data.k2 +  ...
      y1.ter * data.et1 * data.ez2            * data.k2 * data.k2;  
%  magnetic double layer potential
DL2 = cross( data.ez2, y2( 3 ).tez1, k2 ) / z / data.k1 * data.k2 +  ...
      cross( data.ez1, y2( 3 ).tmz2, k1 )                         -  ...
      y1.ter * data.ez1 * data.et2        / z * data.k1 * data.k2 +  ...
      y1.tmr * data.et1 * data.ez2            * data.k2 * data.k2;      
    
%  save output, multiply with prefactors
data = struct;
data.SL1 = double( SL1, [ i1, k1, i2, k2 ] ) * mu2;
data.SL2 = double( SL2, [ i1, k1, i2, k2 ] ) * eps2;
data.DL1 = double( DL1, [ i1, k1, i2, k2 ] );
data.DL2 = double( DL2, [ i1, k1, i2, k2 ] );
