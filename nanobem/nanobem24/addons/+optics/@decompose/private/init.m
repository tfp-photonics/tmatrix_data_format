function obj = init( obj, k0, varargin )
%  Initialize plane wave decomposition.
%
%  Usage for obj = optics.decompose :
%    obj = init( obj, k0, pol, dir )
%  Input
%    k0       :  wavenumber of light in vacuum
%    efield   :  electric field components
%    dir      :  propagation directions

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'efield', [] );
addOptional( p, 'dir', [] );
%  parse input
parse( p, varargin{ : } );

obj.k0 = k0;
obj.efield = p.Results.efield;
obj.dir = p.Results.dir;
