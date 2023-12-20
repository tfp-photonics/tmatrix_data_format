function obj = init( obj, layer, pt, varargin )
%  INIT - Initialize dipole object.
%
%  Usage for obj = stratified.dipole :
%    obj = init( layer, layer, pt, varargin )
%  Input
%    layer  :  layer structure
%    pt     :  dipole positions

obj.layer = layer;
obj.pt = pt;
%  galerkin.dipole object and options 
obj.dip = galerkin.dipole( pt, varargin{ : } );
obj.opt = varargin;  
