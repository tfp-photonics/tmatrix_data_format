function obj = convert( obj, key )
%  CONVERT - Convert between different definitions of T-matrices.
%
%  Usage for obj = multipole.tmatrix :
%    obj = convert( obj, key )
%  Input
%    key    :  'to_h5' or 'from_h5'

for i = 1 : numel( obj )
  switch key
    case 'to_h5'
      obj( i ).aa =        obj( i ).aa;
      obj( i ).ab = - 1i * obj( i ).ab;
      obj( i ).ba = + 1i * obj( i ).ba;
      obj( i ).bb =        obj( i ).bb;      
    case 'from_h5'
      obj( i ).aa =         obj( i ).aa;
      obj( i ).ab =  + 1i * obj( i ).ab;
      obj( i ).ba =  - 1i * obj( i ).ba;
      obj( i ).bb =         obj( i ).bb;         
  end
end
