function obj = init( obj, solver, varargin )
%  INIT - Initialize T-matrix.
%
%  Usage for obj = multipole.tmatrix :
%    obj = init( obj, solver, k0 )
%  Input
%    solver   :  T-matrix or Mie solver
%    tau      :  boundary elements
%    k0       :  wavenumber of light in vacuum
      
%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'k0', [] );
%  parse input
parse( p, varargin{ : } );

obj.solver = solver;
obj.k0 = p.Results.k0;
