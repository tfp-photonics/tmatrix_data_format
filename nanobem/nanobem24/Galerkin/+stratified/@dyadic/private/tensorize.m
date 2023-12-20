function [ y1, y2, data ] = tensorize( tab, pos1, pos2, k0, varargin )
%  TENSORIZE - Green functions and auxiliary data as tensor objects.
%
%  Usage :
%    [ y1, y2, data ] = tensorize( tab, pos1, pos2, k0, PropertyPairs )
%  Input
%    tab    :  tabulated Green function object
%    pos1   :  field positions
%    pos2   :  source positions
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    ind    :  tensor indices [i1,i2,k1,k2]
%  Output
%    y1     :  Green function values, tensor
%    y2     :  derivatives of Green functions, tensor
%    data   :  auxiliary material properties

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'ind', 1 : 4 );
%  parse input
parse( p, varargin{ : } );

%  Green functions and derivatives
[ y1, y2 ] = interp2( tab, pos1, pos2, k0, varargin{ : } );
%  dummy indices for tensor class
ind = num2cell( p.Results.ind );
[ i1, i2, k1, k2 ] = deal( ind{ : } );
%  convert to tensors
for name = convertCharsToStrings( fieldnames( y1 ) ) .'
  %  Green function value
  y1.( name ) = tensor( y1.( name ), [ i1, i2 ] );
  %  derivatives
  y2( 1 ).( name ) = tensor( y2( 1 ).( name ), [ i1, i2, k1 ] );
  y2( 2 ).( name ) = tensor( y2( 2 ).( name ), [ i1, i2, k2 ] );
  y2( 3 ).( name ) = tensor( y2( 3 ).( name ), [ i1, k1, i2, k2 ] );
end

%  radial distance
pos1( :, 3 ) = 0;
pos2( :, 3 ) = 0;
r = pdist2( pos1, pos2 );
%  avoid division by zero
r( r < 1e-10 ) = 1e-10;

%  unit vectors in radial direction
data.r = tensor( r, [ i1, i2 ] );
data.er1 = ( tensor( pos1, [ i1, k1 ] ) - tensor( pos2, [ i2, k1 ] ) ) ./ data.r;
data.er2 = ( tensor( pos1, [ i1, k2 ] ) - tensor( pos2, [ i2, k2 ] ) ) ./ data.r;
%  unit vectors in Cartesian directions
data.ex1 = tensor( [ 1, 0, 0 ], k1 );  data.ex2 = tensor( [ 1, 0, 0 ], k2 ); 
data.ey1 = tensor( [ 0, 1, 0 ], k1 );  data.ey2 = tensor( [ 0, 1, 0 ], k2 );
data.ez1 = tensor( [ 0, 0, 1 ], k1 );  data.ez2 = tensor( [ 0, 0, 1 ], k2 );
%  unit vector in azimuthal direction
data.et1 = cross( data.ez1, data.er1, k1 ); 
data.et2 = cross( data.ez2, data.er2, k2 ); 

%  materials of layer material
ind = indlayer( tab );
mat1 = tab.layer.mat( ind( 1 ) );
mat2 = tab.layer.mat( ind( 2 ) );
%  material properties of layer materials
[ data.k1, data.mu1, data.eps1 ] = deal( mat1.k( k0 ), mat1.mu( k0 ), mat1.eps( k0 ) );
[ data.k2, data.mu2, data.eps2 ] = deal( mat2.k( k0 ), mat2.mu( k0 ), mat2.eps( k0 ) );
