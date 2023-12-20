function pt = Point( layer, tau, pos )
%  POINT - Place points in dielectric environment of layer structure.
%
%  Usage :
%    pt = stratified.Point( layer, tau, pos )
%  Input
%    layer    :  layer structure
%    tau      :  boundary elements
%    pos      :  point positions
%  Output
%    pt       :  points embedded in dielectric environment

%  boundary elements connected to layer structure
inout = vertcat( tau.inout );
ind = inout( :, 2 ) <= layer.n + 1;
%  use same index for all layer materials
for it = find( ind ) .'
  tau( it ).inout = [ inout( it, 1 ), 1 ];
end
%  place points in dielectric environment
pt = Point( tau, pos, 1 : layer.n );

%  material indices of points in layer media
ind = vertcat( pt.imat ) == 1;
imat = indlayer( layer, vertcat( pt.pos ) );
%  set material indices of points
for it = find( ind ) .'
  pt( it ).imat = imat( it );
end
