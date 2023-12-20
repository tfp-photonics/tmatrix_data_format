function [ obj, wait ] = refine1( obj, wait )
%  REFINE1 - Refine quasistatic terms.
%
%  Usage for obj = stratified.pot1.intra2 :
%    [ obj, wait ] = refine1( obj, wait )
%  Input
%    wait   :  waitbar object

%  Green function object and allocate output
green = stratified.tab.intra2( obj.layer, obj.i1, [], [], [] );
yout = cell( numel( obj.yout ), 2 );
%  loop over refinement integrators
for it = 1 : numel( obj.yout )
  %  boundary integrators and quadrature points
  [ oint, data ] = oint2( obj, it );
  y1 = quasistatic( green, data.pos1, data.pos2, [], 'rows', 1 );
  
  %  loop over interfaces
  for i1 = 1 : 2
    %  multiply SL terms with shape functions
    yy.zp = oint.zp( y1( i1 ).z  );
    yy.pz = oint.pz( y1( i1 ).z  );
    yy.pp = oint.pp( y1( i1 ).zz );
    yy.zz = oint.zz( y1( i1 ).zz );
    yy.ss = oint.ss( y1( i1 ).s  );
    %  multiply DL terms with shape functions
    yy.rp = oint.rp( y1( i1 ).rz );
    yy.pr = oint.pr( y1( i1 ).rz );  
    yy.rz = oint.rz( y1( i1 ).r  );
    yy.zr = oint.zr( y1( i1 ).r  );  
  
    %  indices
    yy.i1 = data.i1;
    yy.i2 = data.i2;
    %  save refined elements
    yout{ it, i1 } = yy;
  end
  
  %  update Green function object
  wait = stratified.pot1.waitbar( 'update', wait );
end

%  vectorize output
yout1 = horzcat( yout{ :, 1 } );
yout2 = horzcat( yout{ :, 2 } );
%  concatenate refined elements
obj.yout = struct;
for name = convertCharsToStrings( fieldnames( yout1 ) ) .'
  obj.yout( 1 ).( name ) = cat( 1, yout1.( name ) );
  obj.yout( 2 ).( name ) = cat( 1, yout2.( name ) );
end
