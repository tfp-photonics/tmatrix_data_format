function pts = set( obj, tau1, tau2, varargin )
%  SET - Set boundary elements and precomute SL and DL elements.
%
%  Usage for obj = galerkin.pot1.base :
%    obj = set( obj, tau1, tau2, PropertyPairs )
%  Input
%    tau1   :  first set of boundary elements
%    tau2   :  second set of boundary elements
%  PropertyName
%    'waitbar'  :  use waitbar during initialization
%    'memax'    :  slice integrators to bunches of maximal size MEMAX

%  slice boundary elements into bunches of equal shape and material
tau1 = slice( tau1 );
tau2 = slice( tau2 );
%  allocate array for integration points
pts = cell( numel( tau1 ), numel( tau2 ) );

%  assert that first integrator is of type STD
assert( isa( obj( 1 ), 'galerkin.pot1.std' ) );
%  initialzize STD integrators
for i1 = 1 : numel( tau1 )
for i2 = 1 : numel( tau2 )
  if ~isempty( intersect( tau1{ i1 }( 1 ).inout, tau2{ i2 }( 1 ).inout ) )
    pts{ i1, i2 } = set1( obj( 1 ), tau1{ i1 }, tau2{ i2 } );
  end
end
end
%  slice into bunches of maximal size MEMAX
pts = num2cell( slice( horzcat( pts{ : } ), varargin{ : } ) );

%  initialize refinement integrators
for i1 = 1 : numel( pts )
  %  STD integrator
  pt1 = pts{ i1 };
  %  combine with refinement integrators
  pt2 = arrayfun(  ...
    @( x ) set1( x, pt1.tau1, pt1.tau2 ), obj( 2 : end ), 'uniform', 0 );
  pts{ i1 } = [ pt1, horzcat( pt2{ : } ) ];
end
%  slice integrator array into bunches of maximal size MEMAX
pts = slice( horzcat( pts{ : } ), varargin{ : } );

%  precompute elements for evaluation of SL and DL potential
pts = feval1( pts, 'eval1',  ...
  horzcat( pts.initialize ), 'name', 'galerkin.pot1.base', varargin{ : } );
