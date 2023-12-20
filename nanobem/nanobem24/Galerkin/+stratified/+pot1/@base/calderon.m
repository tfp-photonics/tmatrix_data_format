function cal = calderon( obj, green, k0, varargin )
%  CALDERON - Compute reflected Calderon matrix.
%
%  Usage for obj = stratified.pot1.base :
%    cal = calderon( obj, green, k0, PropertyPairs )
%  Input
%    k0     :  wavenumber of light in vacuum
%    green  :  tabulated Green function object
%  Output
%   cal     :  reflected Calderon matrix

%  compute reflected SL and DL potentials
data = eval( obj, green, k0, varargin{ : } );
%  reflected Calderon submatrices
A11 = data.DL1;
A22 = data.DL2;
A12 = - 1i * k0 * data.SL1;
A21 =   1i * k0 * data.SL2;
%  assemble Calderon matrix
cal = [ A11, A12; A21, A22 ];
