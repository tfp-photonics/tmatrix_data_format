function obj = init( obj, mat, z )
%  INIT - Initialize layer structure.
%
%  Usage for obj = stratified.layerstructure :
%    obj = init( obj, mat, z )
%  Input
%    mat  :  material properties
%    z    :  z-values of interfaces

%  save input
[ obj.mat, obj.z ] = deal( mat, z );
%  layer thickness
obj.d = z( 2 : end ) - z( 1 : end - 1 );
