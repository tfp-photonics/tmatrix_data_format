function cal = calderon( obj, k0 )
%  CALDERON - Calderon matrix for contour integral method.
%
%  Usage for obj = feibelman.bemsolver :
%    cal = calderon( obj, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum
%  Output
%    cal    :  Calderon matrix

%  boundary matrices
[ b1, b2, I ] = bmatrix( obj, k0 );
%  compute SL and DL potential at particle inside and outside
data1 = eval( obj.pot, k0, 1 );
data2 = eval( obj.pot, k0, 2 );
%  matrices for Calderon identities
id1 = 0.5 * I - [ data1.DL, - 1i * k0 * data1.SL1; 1i * k0 * data1.SL2, data1.DL ];
id2 = 0.5 * I + [ data2.DL, - 1i * k0 * data2.SL1; 1i * k0 * data2.SL2, data2.DL ];

%  relate U1 to U2 using boundary matrices
for it = 1 : numel( obj.ind )
  i1 = [ obj.ind{ it }; ndof( obj.tau ) + obj.ind{ it } ];
  id1( i1, i1 ) = ( id1( i1, i1 ) / full( b1( i1, i1 ) ) ) * b2( i1, i1 );
end
%  Calderon matrix
cal = id2 - id1;
