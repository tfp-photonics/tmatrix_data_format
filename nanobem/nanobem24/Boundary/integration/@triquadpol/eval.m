function [ pos, w, f, fp ] = eval( obj )
%  EVAL - Evaluate integration points and basis functions.
%
%  Usage for obj = triquadpol :
%    [ pos, w, f, fp ] = eval( obj )
%  Output
%    pos    :  integration positions
%    w      :  integration weights
%    f      :  shape functions
%    fp     :  divergence of shape functions (only for edge elements)

%  extract boundary elements and quadrature rules
[ tau, quad ] = deal( obj.tau( obj.i1 ), obj.quad );
%  evaluate basis functions and quadrature points
[ f, fp, pos ] = basis( tau, quad.x, quad.y, 'same', 1 );
w = quad.w;
