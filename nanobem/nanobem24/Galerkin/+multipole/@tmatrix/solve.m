function sol = solve( obj, q )
%  SOLVE - Solve T-matrix equation for incoming multipole coefficients.
%
%  Usage for obj = multipole.tmatrix :
%    sol = solve( obj, q )
%  Input
%    q    :  incoming multipole coefficients
%  Output
%    sol  :  multipole solution

%  multipole coefficients
a = obj.aa * q.a + obj.ab * q.b;
b = obj.ba * q.a + obj.bb * q.b;
%  set output
solver = obj.solver;
sol = multipole.solution(  ...
  solver.tab, solver.embedding, obj.k0, a, b, q.a, q.b );
