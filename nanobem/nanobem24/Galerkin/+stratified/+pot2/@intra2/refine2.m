function [ y1, y2, data ] = refine2( obj, data1, k0, varargin )
%  REFINE2 - Handle refined quasistatic terms.
%
%  Usage for obj = stratified.pot2.intra2 :
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
green = stratified.tab.intra2( obj.layer, obj.i1, [], [], [] );
[ y1, y2 ] = interp1( obj, green,  ...
   data1.pos1, data1.pos2, k0, 'stat', 1, 'refl', 1, varargin{ : } );  

%  deal with refined Green function elements
if isempty( obj.yout ) || isempty( y1 )
  data = [];
else
  %  set refined Green function elements to zero
  yout = obj.yout;
  [ y1, siz ] = fzero( obj, y1, p.Results.ind( 1 : 3 ), yout( 1 ).i1, yout( 1 ).i2 );
  [ y2, ~   ] = fzero( obj, y2, p.Results.ind( 1 : 4 ), yout( 1 ).i1, yout( 1 ).i2 );
  
  %  quasistatic reflection coefficients
  rq = fresnelquasistatic( green, k0 );
  %  allocate single and double layer potentials
  [ SL1, SL2, DL1, DL2 ] = deal( 0 );  
  
  %  loop over interfaces
  for m = 1 : 2  
    %  refined electric single layer elements
    SL1 = SL1 - ( rq( m ).tm + rq( m ).te ) * yout( m ).x  -  ...
                  rq( m ).tm / data1.k1 ^ 2 * yout( m ).pp -  ...
                  rq( m ).te * yout( m ).zz + rq( m ).te * yout( m ).ss;
    %  refined magnetic single layer elements
    SL2 = SL2 - ( rq( m ).tm + rq( m ).te ) * yout( m ).x  -  ...
                  rq( m ).te / data1.k1 ^ 2 * yout( m ).pp -  ...
                  rq( m ).tm * yout( m ).zz + rq( m ).te * yout( m ).ss;
           
    %  refined electric double layer potentials
    DL1 = DL1 +  ...
      rq( m ).tm * yout( m ).pr - rq( m ).tm * yout( m ).zr * data1.k1 ^ 2 -  ...
      rq( m ).te * yout( m ).rp + rq( m ).te * yout( m ).rz * data1.k1 ^ 2;
    %  refined magnetic double layer potentials
    DL2 = DL2 +  ...
      rq( m ).te * yout( m ).pr - rq( m ).te * yout( m ).zr * data1.k1 ^ 2 -  ...
      rq( m ).tm * yout( m ).rp + rq( m ).tm * yout( m ).rz * data1.k1 ^ 2;      
  end
  
  %  allocate output
  [ data.SL1, data.SL2, data.DL1, data.DL2 ] =  ...
                   deal( zeros( siz( 1 ) * siz( 2 ), size( SL1, 2 ), 3 ) );
  %  set refined elements
  ind = sub2ind( siz( [ 1, 2 ] ), yout( 1 ).i1, yout( 1 ).i2 );
  data.SL1( ind, :, : ) = SL1;  
  data.SL2( ind, :, : ) = SL2;  
  data.DL1( ind, :, : ) = DL1;  
  data.DL2( ind, :, : ) = DL2;  
  %  reshape output
  for name = [ "SL1", "SL2", "DL1", "DL2" ]
    data.( name ) = reshape( data.( name ), siz( 1 ), siz( 2 ), [], 3 );
  end
end
