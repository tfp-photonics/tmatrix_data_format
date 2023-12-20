function [ d1, d2 ] = eval2( obj, tau, k0 )
%  EVAL - Evaluate Feibelman parameters.
%
%  Usage for obj = feibelman.param :
%    obj = eval1( obj, tau, k0 )
%  Input
%    tau    :  boundary elements
%    k0     :  wavenumber of light in vacuum

%  assign output
[ d1, d2 ] = deal( zeros( numel( tau ), 1 ) );
%  centroid positions of boundary elements
pos = vertcat( tau.pos );

%  loop over Feibelman parameters
for i = 1 : numel( obj )
  %  evaluate object
  [ e1, e2 ] = eval1( obj( i ), k0 );
  %  assign parameters
  [ i1, ind ] = ismember( vertcat( obj( i ).tau.pos ), pos, 'rows' );
  d1( ind( i1 ) ) = e1( i1 );
  d2( ind( i1 ) ) = e2( i1 );
end
