function sol = match( obj, sol, varargin )
%  MATCH - Match tangential fields from boundary outside to inside.
%
%  Usage for obj = feibelman.bemsolver :
%    sol = solve( obj, sol )
%    sol = solve( obj, sol, b1, b2 )
%  Input
%    sol    :  galerkin.solution object with fields at particle outside
%    b1     :  boundary matrix at particle inside
%    b2     :  boundary matrix at particle outside
%  Output
%    sol    :  feibelman.solution object

%  wavenumber of light in vacuum
k0 = sol.k0;
%  boundary matrices
if isempty( varargin )
  [ b1, b2 ] = bmatrix( obj, k0 );
else
  [ b1, b2 ] = deal( varargin{ : } );
end

%  tangential fields at particle outside
n = ndof( obj.tau );
u1 = vertcat( reshape( sol.e, n, [] ), reshape( sol.h, n, [] ) );
%  match fields at particle inside using boundary matrices
for it = 1 : numel( obj.ind )
  i1 = [ obj.ind{ it }; n + obj.ind{ it } ];
  u1( i1, : ) = b1( i1, i1 ) \ ( b2( i1, i1 ) * u1( i1 ) );
end

% set output
sol = feibelman.solution( sol.tau, k0, sol.e, sol.h );
sol.et = reshape( u1( 1 : n,       : ), size( sol.e ) );
sol.ht = reshape( u1( n + 1 : end, : ), size( sol.h ) );
