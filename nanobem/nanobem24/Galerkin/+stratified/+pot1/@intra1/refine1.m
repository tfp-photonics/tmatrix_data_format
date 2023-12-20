function [ obj, wait ] = refine1( obj, wait )
%  REFINE1 - Refine quasistatic terms.
%
%  Usage for obj = stratified.pot1.intra1 :
%    [ obj, wait ] = refine1( obj, wait )
%  Input
%    wait   :  waitbar object

%  Green function object and allocate output
green = stratified.tab.intra1( obj.layer, obj.i1, [], [] );
yout = cell( size( obj.yout ) );
%  loop over refinement integrators
for it = 1 : numel( obj.yout )
  %  boundary integrators and quadrature points
  [ oint, data ] = oint2( obj, it );
  y1 = quasistatic( green, data.pos1, data.pos2, [], 'rows', 1 );
  
  %  multiply SL terms with shape functions
  yy.zp = oint.zp( y1.z  );
  yy.pz = oint.pz( y1.z  );
  yy.pp = oint.pp( y1.zz );
  yy.zz = oint.zz( y1.zz );
  yy.ss = oint.ss( y1.s  );
  %  multiply DL terms with shape functions
  yy.rp = oint.rp( y1.rz );
  yy.pr = oint.pr( y1.rz );  
  yy.rz = oint.rz( y1.r  );
  yy.zr = oint.zr( y1.r  );  
  
  %  indices
  yy.i1 = data.i1;
  yy.i2 = data.i2;
  %  save refined elements
  yout{ it } = yy;
  
  %  update Green function object
  wait = stratified.pot1.waitbar( 'update', wait );
end

%  vectorize output
yout = horzcat( yout{ : } );
%  concatenate refined elements
obj.yout = struct;
for name = convertCharsToStrings( fieldnames( yout ) ) .'
  obj.yout.( name ) = cat( 1, yout.( name ) );
end
