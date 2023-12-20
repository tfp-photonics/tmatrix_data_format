function pot = reflected( layer, tau1, tau2, varargin )
%  REFLECTED - Initialize reflected potential integrator.
%
%  Usage :
%    pot = stratified.pot1.reflected( layer, tau1, tau2, PropertyPairs )
%  Input
%    layer    :  layer structure
%    tau1     :  first set of boundary elements
%    tau2     :  second set of boundary elements
%  PropertyPairs
%    waitbar  :  show waitbar during initialization
%  Ouput
%    pot      :  reflected potential evaluators

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'waitbar', 0 );
addParameter( p, 'name', 'stratified.pot1.reflected' );
%  parse input
parse( p, varargin{ : } );

%  material indices at outside
i1 = vertcat( tau1.inout );  i1 = i1( :, 2 );
i2 = vertcat( tau2.inout );  i2 = i2( :, 2 );
%  allocate output
pot = cell( layer.n + 1);

%  loop over unique material combinations
for m1 = unique( i1( i1 <= layer.n + 1 ) ) .'
for m2 = unique( i2( i2 <= layer.n + 1 ) ) .'
  %  input for potential integrator
  data = [ { layer, tau1( i1 == m1 ), tau2( i2 == m2 ) }, varargin ];
  %  initialize potential integrator
  if m1 == m2 && any( ismember( m1, [ 1, layer.n + 1 ] ) )
    pot{ m1, m2 } = stratified.pot1.intra1( data{ : } );
  elseif m1 == m2
    pot{ m1, m2 } = stratified.pot1.intra2( data{ : } );
  else
    pot{ m1, m2 } = stratified.pot1.inter(  data{ : } );
  end
end
end
%  output vector
pot = horzcat( pot{ : } );

%  initialize waitbar
wait = stratified.pot1.waitbar( 'init', pot, p.Results.waitbar, p.Results.name );
%  loop over potential integrators with refinement
ind = arrayfun( @( x ) ~isempty( x.yout ), pot, 'uniform', 1 );
for it = find( ind )
  [ pot( it ), wait ] = refine1( pot( it ), wait );
end
%  close waitbar
stratified.pot1.waitbar( 'close', wait );
