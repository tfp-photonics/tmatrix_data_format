function [ y, data ] = refine2( obj, data1, k0 )
%  REFINE2 - Handle refined quasistatic terms.
%
%  Usage for obj = stratified.pot1.intra1 :
%    [ y, data ] = refine2( obj, data1, k0 )
%  Input
%    data1  :  positions and material properties
%    k0     :  wavenumber of light in vacuum
%  Output
%    y      :  quasistatic Green function elements
%    data   :  refined single and double layer potential terms

%  quasistatic Green function objects
green = stratified.tab.intra1( obj.layer, obj.i1, [], [] );
y = quasistatic( green, data1.pos1, data1.pos2, k0, 'refl', 1 );

%  deal with refined Green function elements
if isempty( obj.yout ) || isempty( y )
  data = [];
else
  %  set refined Green function elements to zero
  yout = obj.yout;
  names = convertCharsToStrings( fieldnames( y ) ) .';
  for name = names
    siz = size( y.( name ) );
    y1 = reshape(  ...
      permute( y.( name ), [ 1, 3, 2, 4 ] ), [], siz( 2 ), siz( 4 ) );
    y1( sub2ind( siz( [ 1, 3 ] ), yout.i1, yout.i2 ), :, : ) = 0;
    y.( name ) =  ...
      ipermute( reshape( y1, siz( [ 1, 3, 2, 4 ] ) ), [ 1, 3, 2, 4 ] );
  end
  
  %  quasistatic reflection coefficients
  rq = fresnelquasistatic( green, k0 ); 
  %  refined electric single layer elements
  SL1 = rq.tm * yout.pp / data1.k1 ^ 2 -  ...
      ( rq.tm + rq.te ) * ( yout.zp + yout.pz ) -  ...
        rq.te * yout.zz + rq.te * yout.ss;
  %  refined magnetic single layer elements, use duality
  SL2 = rq.te * yout.pp / data1.k1 ^ 2 -  ...
      ( rq.te + rq.tm ) * ( yout.zp + yout.pz ) -  ...
        rq.tm * yout.zz + rq.tm * yout.ss;
      
  %  refined electric double layer potentials
  DL2 = rq.te * ( yout.pr - yout.zr * data1.k1 ^ 2 ) -  ...
        rq.tm * ( yout.rp - yout.rz * data1.k1 ^ 2 );
  %  refined magnetic double layer potentials, use duality
  DL1 = rq.tm * ( yout.pr - yout.zr * data1.k1 ^ 2 ) -  ...
        rq.te * ( yout.rp - yout.rz * data1.k1 ^ 2 );
  
  %  allocate output
  [ data.SL1, data.SL2, data.DL1, data.DL2 ] =  ...
       deal( zeros( siz( 1 ) * siz( 3 ), siz( 2 ), siz( 4 ) ) );
  %  set refined elements
  ind = sub2ind( siz( [ 1, 3 ] ), yout.i1, yout.i2 );
  data.SL1( ind, :, : ) = SL1;  
  data.SL2( ind, :, : ) = SL2;  
  data.DL1( ind, :, : ) = DL1;  
  data.DL2( ind, :, : ) = DL2;  
  %  reshape output
  for name = [ "SL1", "SL2", "DL1", "DL2" ]
    data.( name ) = reshape( data.( name ), siz( [ 1, 3, 2, 4 ] ) );
    data.( name ) = permute( data.( name ), [ 1, 3, 2, 4 ] );
  end
end

