function data = eval( obj, green, k0, varargin )
%  EVAL - Evaluate reflected single and double layer potentials.
%
%  Usage for obj = stratified.pot1.base :
%    data = eval( obj, green, k0, PropertyPairs )
%  Input
%    green    :  tabulated Green function object
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    siz      :  size of output matrices
%  Output
%    data.SL  :  reflected single layer potential
%    data.DL  :  reflected double layer potential

%  maxiumum number of evaluation points and degrees of freedom
n1 = max( arrayfun( @( x ) max( x.pt1.nu ), obj, 'uniform', 1 ) );
n2 = max( arrayfun( @( x ) ndof( x.tau2 ),  obj, 'uniform', 1 ) );
%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'siz', [ n1, 3, n2 ] );
parse( p, varargin{ : } );

%  allocate output
[ data.SL1, data.SL2, data.DL1, data.DL2 ] = deal( 0 );
%  loop over reflected potential evaluators
for it = 1 : numel( obj )
  %  evaluate reflected potential
  data1 = eval1( obj( it ), green, k0 );
  
  %  conversion between local and global degrees of freedom
  nu1 = obj( it ).pt1.nu( : );
  nu2 = vertcat( obj( it ).tau2.nu );
    %  accumulate matrix
  [ nu1, nu2, k ] = ndgrid( nu1( : ), nu2( : ), 1 : 3 );
  fun = @( x ) accumarray( {  ...
    nu1( : ), nu2( : ), k( : ) }, x( : ), p.Results.siz( [ 1, 3, 2 ] ) );
  %  update single and double layer potentials
  %    DL potential has different signs in Chew and Hohenester, Eq. (11.39)
  data.SL1 = data.SL1 + fun( data1.SL1 );
  data.SL2 = data.SL2 + fun( data1.SL2 );
  data.DL1 = data.DL1 - fun( data1.DL1 );
  data.DL2 = data.DL2 - fun( data1.DL2 );
end

%  reshape output
for name = convertCharsToStrings( fieldnames( data ) ) .'
  data.( name ) = permute( data.( name ), [ 1, 3, 2 ] );
end
