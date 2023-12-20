function data = eval( obj, pos1, pos2, k0, varargin )
%  EVAL - Evaluate dyadic Green functions.
%
%  Usage for obj = stratified.dyadic : 
%    data = eval( obj, pos1, pos2, k0, PropertyPairs );
%  Input
%    tab      :  tabulated Green function
%    pos1     :  observation points
%    pos2     :  source points
%    k0       :  wavenumber of light in vacuum
%  Output
%    data     :  single and double layer potentials

%  layer indices for positions
ind1 = indlayer( obj.tab( 1 ).layer, pos1 );
ind2 = indlayer( obj.tab( 1 ).layer, pos2 );
%  layer indices of Green function objects
ind = arrayfun( @( x ) indlayer( x ), obj.tab, 'uniform', 0 );
ind = vertcat( ind{ : } );

%  allocate output
[ data.SL1, data.DL1, data.SL2, data.DL2 ] =  ...
               deal( zeros( [ size( pos1 ), size( pos2 ) ] ) );
             
%  loop over unique indices   
for i1 = unique( reshape( ind1, 1, [] ) )
for i2 = unique( reshape( ind2, 1, [] ) )
  %  find Green function object for layer indices
  [ ~, it ] = ismember( [ i1, i2 ], ind, 'rows' );
  %  evaluate Green function  
  L1 = ind1 == i1;
  L2 = ind2 == i2;
  if i1 == i2
    data1 = eval1( obj.tab( it ), pos1( L1, : ), pos2( L2, : ), k0, varargin{ : } );
  else
    data1 = eval2( obj.tab( it ), pos1( L1, : ), pos2( L2, : ), k0, varargin{ : } );
  end
  %  set output
  for name = [ "SL1", "DL1", "SL2", "DL2" ]
    data.( name )( L1, :, L2, : ) = data1.( name );
  end
end
end

