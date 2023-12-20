function obj = init( obj, mat, imat, varargin )
%  INIT - Initialize T-matrix solver.
%
%  Usage for obj = multipole.tsolver :
%    obj = init( obj, mat, imat, PropertyPairs )
%  Input
%    mat      :  material vector
%    imat     :  index for embedding medium      
%  PropertyName
%    rules    :  quadrature rules

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules', quadboundary.rules );
%  parse input
parse( p, varargin{ : } );

obj.mat  = mat;
obj.imat = imat;
obj.rules = p.Results.rules;
