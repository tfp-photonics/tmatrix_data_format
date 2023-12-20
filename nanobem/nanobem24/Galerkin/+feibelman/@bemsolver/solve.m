function [ sol1, sol2 ] = solve( obj, qinc )
%  SOLVE - Solve BEM equations.
%
%  Usage for obj = feibelman.bemsolver :
%    sol = solve( obj, qinc )
%  Input
%    qinc   :  structure with inhomogeneities and wavenumber
%  Output
%    sol1   :  solution with Feibelman parameters
%    sol2   :  solution w/o  Feibelman parameters

%  number of global degrees of freedom
tau = obj.tau;
n = ndof( tau );
%  boundary matrices
k0 = qinc.k0;
[ b1, b2, I ] = bmatrix( obj, k0 );

%  compute SL and DL potential at particle inside and outside
data1 = eval( obj.pot, k0, 1 );
data2 = eval( obj.pot, k0, 2 );
%  matrices for Calderon identities
id1 = 0.5 * I - [ data1.DL, - 1i * k0 * data1.SL1; 1i * k0 * data1.SL2, data1.DL ];
id2 = 0.5 * I + [ data2.DL, - 1i * k0 * data2.SL1; 1i * k0 * data2.SL2, data2.DL ];

%  inhomogeneity for incoming fields 
q = vertcat( reshape( qinc.e, n, [] ), reshape( qinc.h, n, [] ) );
%  solution w/o Feibelman parameters (only if requested)
if nargout == 2
  u = ( id2 - id1 ) \ q;
  sol2 = galerkin.solution( tau, k0, u( 1 : n, : ), u( n + 1 : end, : ) );
  sol2.e = reshape( sol2.e, size( qinc.e ) );
  sol2.h = reshape( sol2.h, size( qinc.h ) );
end

%  relate U1 to U2 using boundary matrices
for it = 1 : numel( obj.ind )
  i1 = [ obj.ind{ it }; n + obj.ind{ it } ];
  id1( i1, i1 ) = ( id1( i1, i1 ) / full( b1( i1, i1 ) ) ) * b2( i1, i1 );
end
%  solution with Feibelman parameters
u2 = ( id2 - id1 ) \ q;
sol1 = galerkin.solution( tau, k0, u2( 1 : n, : ), u2( n + 1 : end, : ) );
sol1.e = reshape( sol1.e, size( qinc.e ) );
sol1.h = reshape( sol1.h, size( qinc.h ) );
%  match fields at particle inside using boundary matrices
sol1 = match( obj, sol1, b1, b2 );
