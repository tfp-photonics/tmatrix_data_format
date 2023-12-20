function ind = indlayer( obj, pos )
%  INDLAYER - Locate positions in layer structure.
%
%  Usage for obj = stratified.layerstructure :
%    ind = indlayer( obj, pos )
%  Input
%    pos    :  positions
%  Output
%    ind    :  layer index

%  bin edges
edges = [ - inf, reshape( obj.z, 1, [] ), inf ];
%  layer indices
[ ~, ~, ind ] = histcounts( pos( :, 3 ), edges );
