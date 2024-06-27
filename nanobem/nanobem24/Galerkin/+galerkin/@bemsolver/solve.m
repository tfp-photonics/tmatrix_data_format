function [ sol, obj ] = solve( obj, q )
%  SOLVE - Solve BEM equations.
%
%  Usage for obj = galerkin.bemsolver :
%    sol = solve( obj, q )
%  Input
%    q      :  structure with inhomogeneities and wavenumber
%  Output
%    sol    :  solution with tangential electromagnetic fields

%  inhomogeneity vector
n = ndof( obj.tau );
vec = reshape( cat( 1, q.e, q.h ), 2 * n, [] );
%  solve BEM equations 
obj = fill( obj, q.k0 );
u = obj.cal \ vec;
%  tangential electric and magnetic field
ue = reshape( u(     1 : n,   : ), size( q.e ) );
uh = reshape( u( n + 1 : end, : ), size( q.h ) );
%  set output
sol = galerkin.solution( obj.tau, q.k0, ue, uh );
