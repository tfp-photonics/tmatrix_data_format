function [ y, data ] = refine2( obj, data1, k0 )
%  REFINE2 - Handle refined quasistatic terms.
%
%  Usage for obj = stratified.pot1.inter :
%    [ y, data ] = refine2( obj, data1, k0 )
%  Input
%    data1  :  positions and material properties
%    k0     :  wavenumber of light in vacuum
%  Output
%    y      :  quasistatic Green function elements
%    data   :  refined single and double layer potential terms

%  quasistatic Green function objects
green = stratified.tab.inter( obj.layer, obj.i1, obj.i2, [], [], [] );
y = quasistatic( green, data1.pos1, data1.pos2, k0, 'refl', 1 );

%  deal with refined Green function elements
if isempty( obj.yout )
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

  %  quasistatic transmission coefficients
  tq = fresnelquasistatic( green, k0 );
  %  wavenumbers 
  k1 = data1.k1;
  k2 = data1.k2;
  %  ratio of impedances
  z = sqrt( data1.mu1 / data1.eps1 ) / sqrt( data1.mu2 / data1.eps2 );
    
  %  refined electric single layer elements
  SL1 = tq.tm * z *   yout.pp  / k1 / k2 -                         ...
        tq.tm * z * ( yout.zp2 * k1 / k2 + yout.pz1 / k1 * k2 ) -  ...
        tq.te     * ( yout.zp1 + yout.pz2 ) -                      ...
        tq.te * yout.zz + tq.te * yout.ss;
  %  refined magnetic single layer elements
  SL2 = tq.te / z *   yout.pp  / k1 / k2    -                      ...
        tq.te / z * ( yout.zp2 * k1 / k2 + yout.pz1 / k1 * k2 ) -  ...
        tq.tm     * ( yout.zp1 + yout.pz2 ) -                      ...
        tq.tm * yout.zz + tq.tm * yout.ss;

  %  refined magnetic double layer potentials
  DL1 = tq.tm * z * ( yout.pr * k2 / k1 - yout.zr * k1 * k2 ) -  ...
        tq.te     * ( yout.rp           - yout.rz * k2 * k2 );      
  %  refined double layer potentials
  DL2 = tq.te / z * ( yout.pr * k2 / k1 - yout.zr * k1 * k2 ) -  ...
        tq.tm     * ( yout.rp           - yout.rz * k2 * k2 );
  
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

