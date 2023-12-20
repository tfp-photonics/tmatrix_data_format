function [ y, siz ] = fzero( ~, y, ind, i1, i2 )
%  FZERO - Set indexed quasistatic terms to zero.
%
%  Usage for obj = stratified.pot2.base :
%    [ y, siz ] = fzero( obj, y, ind, i1, i2 )
%  Input
%    y    :  quasistatic terms
%    ind  :  indices for tensor class
%    i1   :  index for observation points
%    i2   :  index for source points
%  Output
%    y    :  updated quasistatic terms 
%    siz  :  size of arrays

for name = convertCharsToStrings( fieldnames( y ) ) .'
  %  get elements
  yy = double( y.( name ), ind );
  %  reshape input
  siz = size( yy );
  yy = reshape( yy, siz( 1 ) * siz( 2 ), [] );
  %  set indexed elements to zero
  yy( sub2ind( siz( [ 1, 2 ] ), i1, i2 ), : ) = 0;
  yy = reshape( yy, siz );  
  %  convert output to tensor class
  y.( name ) = tensor( yy, ind );
end
