function [ q, data ] = eval( obj, tau, k0, varargin )
%  EVAL - Inhomogeneities for dipole excitation.
%
%  Usage for obj = stratified.dipole :
%    [ q, data ] = eval( obj, tau, k0 )
%  Input
%    tau    :  boundary elements
%    k0     :  wavelength of light in vacuum
%  Output
%    q      :  structure containing inhomogeneities for BEM
%    data   :  potential integrators and Green function for reuse

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'pot',   [] );
addParameter( p, 'green', [] );
%  parse input
parse( p, varargin{ : } );

%  extract input
[ pot, green ] = deal( p.Results.pot, p.Results.green );
%  reflected potential integrator
if isempty( pot )
  pot = stratified.pot2.reflected( obj.layer, obj.pt, tau, obj.opt{ : } );
end
%  reflected Green function
if isempty( green )
  pos1 = reshape( eval( quadboundary( tau ) ), [], 3 );
  r = slice( obj.layer, pos1, vertcat( obj.pt.pos ) );
  r = grid( obj.layer, r, obj.opt{ : } );
  green = stratified.tab.green( r, obj.opt{ : } );
end

%  fill Green function
green = fill( green, k0 );
%  number of points and degrees of freedom 
n1 = numel( obj.pt );
n2 = ndof( tau );
%  allocate output
e = zeros( n2, n1, 3 );
h = zeros( n2, n1, 3 );
%  direct excitation
q = eval( obj.dip, tau, k0 );
e( :, 1 : size( q.e, 2 ), : ) = q.e;
h( :, 1 : size( q.h, 2 ), : ) = q.h;

for it = 1 : numel( pot )
  %  evaluate SL and DL potential, evaluation points
  data = eval1( pot( it ), green, k0 );
  [ pt1, pt2 ] = deal( pot( it ).pt1, pot( it ).pt2 );  
    
  %  indices for evaluation points and edges
  [ nu1, nu2 ] = ndgrid( pt1.nu, vertcat( pt2.tau.nu ) );
  subs = { nu1( : ), nu2( : ) }; 
  %  accumulation function
  SL = cat( 3,  ...
    accumarray( subs, reshape( data.SL1( :, :, :, 1 ), [], 1 ), [ n1, n2 ] ),  ...
    accumarray( subs, reshape( data.SL1( :, :, :, 2 ), [], 1 ), [ n1, n2 ] ),  ...
    accumarray( subs, reshape( data.SL1( :, :, :, 3 ), [], 1 ), [ n1, n2 ] ) );
  DL = cat( 3,  ...
    accumarray( subs, reshape( data.DL2( :, :, :, 1 ), [], 1 ), [ n1, n2 ] ),  ...
    accumarray( subs, reshape( data.DL2( :, :, :, 2 ), [], 1 ), [ n1, n2 ] ),  ...
    accumarray( subs, reshape( data.DL2( :, :, :, 3 ), [], 1 ), [ n1, n2 ] ) ); 

  %  inhomogeneities
  e = e +  k0 ^ 2 * permute( SL, [ 2, 1, 3 ] );
  h = h + 1i * k0 * permute( DL, [ 2, 1, 3 ] );
end  

%  set output
q = struct( 'e', e, 'h', h, 'tau', tau, 'k0', k0 );
data = struct( 'pot', pot, 'green', green );
