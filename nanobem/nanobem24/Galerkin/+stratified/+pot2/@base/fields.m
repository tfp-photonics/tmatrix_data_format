function [ e, h ] = fields( obj, green, sol, varargin )
%  FIELDS - Reflected electromagnetic fields at point positions.
%
%  Usage for obj = stratified.pot2.base :
%    [ e, h ] = fields( obj, green, sol, PropertyPairs )
%  Input
%    green      :  tabulated Green function object
%    sol        :  solution of BEM equation
%  PropertyName
%    waitbar    :  show waitbar during evaluation
%    n          :  leading dimension of field arrays
%  Output
%    e          :  electric field at requested points
%    h          :  magnetic field at requested points

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'waitbar', 0 );
addParameter( p, 'n', 0 );
addParameter( p, 'name', 'stratified.pot2.base.fields' );
%  parse input
parse( p, varargin{ : } );

%  tangential fields at boundary
ue = reshape( sol.e, size( sol.e, 1 ), [] );
uh = reshape( sol.h, size( sol.h, 1 ), [] );
%  wavenumber of light in vacuum
k0 = sol.k0;

%  allocate output
n = max( arrayfun( @( x ) max( x.pt1.nu ), obj, 'uniform', 1 ) );
n = max( n, p.Results.n );
[ e, h ] = deal( zeros( n, 3, size( ue, 2 ) ) );

for it = 1 : numel( obj )
  %  evaluate SL and DL potential, evaluation points
  data = eval1( obj( it ), green, k0 );
  [ pt1, pt2 ] = deal( obj( it ).pt1, obj( it ).pt2 );
  
  %  single and double layer potential
  SL1 = permute( reshape( data.SL1, numel( pt1.nu ), [], 3 ), [ 1, 3, 2 ] );
  SL2 = permute( reshape( data.SL2, numel( pt1.nu ), [], 3 ), [ 1, 3, 2 ] );
  DL1 = permute( reshape( data.DL1, numel( pt1.nu ), [], 3 ), [ 1, 3, 2 ] );
  DL2 = permute( reshape( data.DL2, numel( pt1.nu ), [], 3 ), [ 1, 3, 2 ] );
  %  multidimensional multiplication function
  fun = @( x, y ) reshape(  ...
                  reshape( x, [], size( y, 1 ) ) * y, size( x, 1 ), 3, [] );
  %  electric and magnetic fields
  e1 = fun( DL1, ue( vertcat( pt2.tau.nu ), : ) );
  e2 = fun( SL2, ue( vertcat( pt2.tau.nu ), : ) );
  h1 = fun( SL1, uh( vertcat( pt2.tau.nu ), : ) );
  h2 = fun( DL2, uh( vertcat( pt2.tau.nu ), : ) );
         
  e( pt1.nu, :, : ) = e( pt1.nu, :, : ) + 1i * k0 * h1 - e1;
  h( pt1.nu, :, : ) = h( pt1.nu, :, : ) - 1i * k0 * e2 - h2; 
end

%  reshape output
siz = size( sol.e );
e = reshape( e, [ p.Results.n, 3, siz( 2 : end ) ] );
h = reshape( h, [ p.Results.n, 3, siz( 2 : end ) ] );
