function [ sol, obj ] = solve( obj, q )
%  SOLVE - Solve BEM equations.
%
%  Usage for stratified.bemsolver :
%    [ sol, obj ] = solve( obj, q )
%  Input
%    q      :  structure with inhomogeneities and wavenumber
%  Output
%    sol    :  solution with tangential electromagnetic fields

%  fill matrices and compute Calderon matrix
obj = fill( obj, q.k0 );
%  solve BEM equations
n = ndof( obj.tau );
u = obj.cal \ reshape( cat( 1, q.e, q.h ), 2 * n, [] );

%  tangential electric and magnetic field
ue = reshape( u(     1 : n,   : ), size( q.e ) );
uh = reshape( u( n + 1 : end, : ), size( q.h ) );
%  set output
sol = stratified.solution( obj.tau, q.k0, ue, uh );
sol.layer = obj.layer;
