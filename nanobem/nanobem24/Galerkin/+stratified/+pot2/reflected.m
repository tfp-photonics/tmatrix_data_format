function pot = reflected( layer, tau1, tau2, varargin )
%  REFLECTED - Initialize reflected potential integrator.
%
%  Usage :
%    pot = stratified.pot2.reflected( layer, tau1, tau2, PropertyPairs )
%  Input
%    layer    :  layer structure
%    tau1     :  evaluation points
%    tau2     :  boundary elements
%  PropertyPairs
%    waitbar  :  show waitbar during initialization
%  Ouput
%    pot      :  reflected potential evaluators

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'waitbar', 0 );
addParameter( p, 'name', 'stratified.pot2.reflected' );
%  parse input
parse( p, varargin{ : } );

%  slice evaluation points and boundary elements
pt1 = iterpoints( tau1 );
pt2 = quadboundary( tau2, varargin{ : } );
%  keep only points and boundary elements connected to layer structure
inout = vertcat( pt2.inout );
pt1 = pt1( horzcat( pt1.imat ) <= layer.n + 1 );
pt2 = pt2( inout( :, 2 ) <= layer.n + 1 );
%  allocate output
pot = cell( numel( pt1 ), numel( pt2 ) );

%  loop over evaluation points and boundary elements
for i1 = 1 : numel( pt1 )
for i2 = 1 : numel( pt2 )
  %  material indices
  m1 = pt1( i1 ).imat;
  m2 = pt2( i2 ).inout( 2 );
  %  input for potential integrator
  data = [ { layer, pt1( i1 ), pt2( i2 ).tau }, varargin ];
  %  initialize potential integrator
  if m1 == m2 && any( ismember( m1, [ 1, layer.n + 1 ] ) )
    pot{ i1, i2 } = stratified.pot2.intra1( data{ : } );  % intra1
  elseif m1 == m2
    pot{ i1, i2 } = stratified.pot2.intra2( data{ : } );  % intra2
  else
    pot{ i1, i2 } = stratified.pot2.inter( data{ : } );   % inter
  end
end
end
%  output vector
pot = horzcat( pot{ : } );
