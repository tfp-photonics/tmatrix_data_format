function [ oint, data ] = oint2( obj, it )
%  OINT2 - Refined boundary integral evaluators.
%
%  Usage for obj = stratified.pot1.base :
%    [ oint, data ] = oint2( obj, it )
%  Input
%    it     :  quadrature points obj.yout(it)
%  Output
%    oint   :  boundary integral evaluators
%    data   :  positions and material properties


%  Duffy and refined integration are treated on the same footing, 
%  we expand the refined integration points to full size

%  quadrature points
pts = obj.yout( it ).pts;
%  evaluate quadrature positions, weights, and basis functions
switch class( pts )
  case 'quadduffy'
    %  integration points and basis functions
    [ pos1, f1, fp1 ] = eval( pts, 1 );
    [ pos2, f2, fp2 ] = eval( pts, 2 );
    %  integration weights and index of refined boundary elements
    w = pts.w;
    [ ~, data.i1 ] = ismember( vertcat( pts.tau1.pos ), vertcat( obj.tau1.pos ), 'rows' );
    [ ~, data.i2 ] = ismember( vertcat( pts.tau2.pos ), vertcat( obj.tau2.pos ), 'rows' );
      
  case 'quadboundary'
    %  integration points, weights, and basis functions
    [ pos1, w1, f1, fp1 ] = eval( pts( 1 ) );
    [ pos2, w2, f2, fp2 ] = eval( pts( 2 ) );
    %  expand positions
    [ i1, i2 ] = ndgrid( 1 : size( w1, 2 ), 1 : size( w2, 2 ) );
    [ i1, i2 ] = deal( i1( : ), i2( : ) );
    pos1 = pos1( :, i1, : );  
    pos2 = pos2( :, i2, : );
    %  expand basis functions
    f1 = f1( :, i1, :, : );  fp1 = fp1( :, i1, : );
    f2 = f2( :, i2, :, : );  fp2 = fp2( :, i2, : );
    %  integration weights  and index of refined boundary elements
    w = w1( :, i1 ) .* w2( :, i2 );
    [ ~, data.i1 ] = ismember( vertcat( pts( 1 ).tau.pos ), vertcat( obj.tau1.pos ), 'rows' );
    [ ~, data.i2 ] = ismember( vertcat( pts( 2 ).tau.pos ), vertcat( obj.tau2.pos ), 'rows' );    
end
  
%  dummy indices for tensor class
[ i, q, a1, a2 ] = deal( 1, 2, 3, 4 );
%  convert shape functions to tensor class
fp1 = tensor( fp1, [ i, q, a1 ] ); 
fp2 = tensor( fp2, [ i, q, a2 ] );
fx1 = tensor( f1( :, :, :, 1 ), [ i, q, a1 ] );
fx2 = tensor( f2( :, :, :, 1 ), [ i, q, a2 ] );
fy1 = tensor( f1( :, :, :, 2 ), [ i, q, a1 ] );
fy2 = tensor( f2( :, :, :, 2 ), [ i, q, a2 ] );  
fz1 = tensor( f1( :, :, :, 3 ), [ i, q, a1 ] );
fz2 = tensor( f2( :, :, :, 3 ), [ i, q, a2 ] );
%  convert position to tensor
xx = tensor( pos1( :, :, 1 ) - pos2( :, :, 1 ), [ i, q ] );
yy = tensor( pos1( :, :, 2 ) - pos2( :, :, 2 ), [ i, q ] );
%  radius
r = sqrt( xx .^ 2 + yy .^ 2 );
r.val( r.val < 1e-10 ) = 1e-10;
%  cross product of shape function with unit vector
fr1 = ( fx1 * yy - fy1 * xx ) ./ r;
fr2 = ( fx2 * yy - fy2 * xx ) ./ r;

%  integration function
ifun = @( f, g )  ...
  double( sum( f * tensor( g .* w, [ i, q ] ), q ), [ i, a1, a2 ] );
%  boundary integral evaluators for single layer potential
oint.pp = @( g ) ifun( fp1 * fp2, g );
oint.zp = @( g ) ifun( fz1 * fp2, g );
oint.pz = @( g ) ifun( fp1 * fz2, g );
oint.zz = @( g ) ifun( fz1 * fz2, g );
oint.ss = @( g ) ifun( fx1 * fx2 + fy1 * fy2, g );
%  boundary integral evaluators for double layer potential
oint.rp = @( g ) ifun( fr1 * fp2, g );
oint.pr = @( g ) ifun( fp1 * fr2, g );
oint.rz = @( g ) ifun( fr1 * fz2, g );
oint.zr = @( g ) ifun( fz1 * fr2, g );

%  save positions
[ data.pos1, data.pos2 ] = deal( pos1, pos2 );
