function [ tot, P0 ] = decayrate( obj, sol, varargin )
%  DECAYRATE - Total decay rate.
%
%  Usage for obj = stratified.dipole :
%    [ tot, P0 ] = decayrate( obj, sol )
%  Input
%    sol    :  solution of BEM equations
%  Output
%    tot    :  enhancement of total decay rate wrt isolated dipole
%    P0     :  dissipated power of isolated dipole

%  direct contribution
[ tot, ~, P0 ] = decayrate( obj.dip, sol );

%  The following code can be slow for multiple dipoles and should be
%  improved in the future.  Here we compute the electric field at the
%  dipole positions for all source dipoles, but only the diagonal
%  components are needed where the source dipole position equals the field
%  dipole position. We could also reuse the tabulated Green functions.

%  reflected dipole and nanoparticle field
e1 = fields( obj, obj.pt, sol.k0, varargin{ : }, 'refl', 1 );
e2 = fields( sol, obj.pt, varargin{ : }, 'refl', 1 );
%  total field
e = e1 + e2;

%  add reflected contributions to decay rate
for i = 1 : numel( obj.pt )
for k = 1 : 3
  tot( i, k ) = tot( i, k ) + 0.5 * sol.k0 * imag( e( i, k, i, k ) ) / P0( i );
end
end
