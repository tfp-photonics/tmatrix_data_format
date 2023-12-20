function obj = init( obj, layer, i1, varargin )
%  INIT - Initialize Sommerfeld integrator.
%
%  Usage for obj = stratified.isommerfeld :
%    obj = init( obj, layer, i1, PropertyPairs )
%  Input
%    layer    :  layer structure
%    i1       :  material index
%  PropertyName
%    semi1    :  scaling factor for real axis of semiellipse
%    semi2    :  scaling factor for imaginary axis of semiellipse 
%    ratio    :  z : r ratio which determines integration path
%    op       :  options for ODE integration      
%    cutoff   :  cutoff parameter for matrix-friendly approach of Chew

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'semi1', 4 );
addParameter( p, 'semi2', 0.05 );
addParameter( p, 'ratio', 2 );
addParameter( p, 'op', odeset( 'AbsTol', 1e-8, 'InitialStep', 1e-3 ) );
addParameter( p, 'cutoff', 0.1 );
%  parse input
parse( p, varargin{ : } );

%  save input
obj.layer  = layer;
obj.i1     = i1;
obj.semi1  = p.Results.semi1;
obj.semi2  = p.Results.semi2;
obj.ratio  = p.Results.ratio;
obj.op     = p.Results.op;
obj.cutoff = p.Results.cutoff;
