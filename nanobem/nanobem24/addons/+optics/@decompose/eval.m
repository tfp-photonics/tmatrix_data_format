function qinc = eval( obj, tau, varargin )
%  EVAL - Inhomogeneity for BEM solver.
%
%  Usage for obj = optics.decompose :
%    qinc = eval( obj, tau, PropertyPairs )
%  Input
%    tau      :  discretized particle boundary
%  PropertyName
%    imat     :  material index for unbounded medium
%    layer    :  layer structure
%    shift    :  additional shift of particle
%    primary  :  reflected fields only
%  Output
%    qinc     :  inhomogeneity for BEM solver

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'shift', [] );
addParameter( p, 'rules', quadboundary.rules );
%  parse input
parse( p, varargin{ : } );

%  shift particle ?
if ~isempty( p.Results.shift )
  for it = 1 : numel( tau )
    tau( it ).verts = tau( it ).verts + p.Results.shift;
    tau( it ).pos   = tau( it ).pos   + p.Results.shift;
  end
end

%  inhomogeneities
[ e, h ] = qbasis( tau,  ...
    @( pt ) fun( obj, pt, varargin{ : } ), p.Results.rules );
%  set output
qinc = struct( 'e', e, 'h', h, 'tau', tau, 'k0', obj.k0 );


function [ e, h ] = fun( obj, pt, varargin )
%  FUN - Electromagnetic fields for planewave decomposition.

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'imat', 1 );
addParameter( p, 'layer', [] );
%  parse input
parse( p, varargin{ : } );

%  positions and material indices
pos = eval( pt );
[ pos, siz ] = deal( reshape( pos, [], 3 ), size( pos ) );
%  allocate output
[ e, h ] = deal( zeros( siz ) );

switch isempty( p.Results.layer )
  case 1    
    %  points connected to medium ?
    if any( pt.inout == p.Results.imat )
      pts = Point( pt.mat, p.Results.imat, pos );
      %  electromagnetic fields
      [ e, h ] = fields1( obj, pts, p.Results.imat );
      [ e, h ] = deal( reshape( e, siz ), reshape( h, siz ) );
    end
  otherwise
    %  points connected to layer structure ?
    layer = p.Results.layer;
    if any( pt.inout <= layer.n + 1 )
      pts = Point( pt.mat, pt.inout( pt.inout <= layer.n + 1 ), pos );
      %  electromagnetic fields
      [ e, h ] = fields2( obj, pts, layer, varargin{ : } );
      [ e, h ] = deal( reshape( e, siz ), reshape( h, siz ) );
    end
end
