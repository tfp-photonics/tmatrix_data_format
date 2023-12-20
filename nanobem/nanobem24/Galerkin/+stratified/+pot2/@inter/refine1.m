function obj = refine1( obj )
%  REFINE1 - Refine quasistatic terms.
%
%  Usage for obj = stratified.pot2.inter :
%    obj = refine1( obj )

%  Green function object and allocate output
green = stratified.tab.inter( obj.layer, obj.i1, obj.i2, [], [], [] );
yout = {};
%  loop over refinement integrators
for it = 1 : numel( obj.yout )
  %  loop over polar integrators
  for pt = obj.yout( it ).pts
    %  boundary integrators and quadrature points
    [ i, q, a, k ] = deal( 1, 2, 3, 4 );
    data = oint2( obj, pt, 'ind', [ i, q, a, k ] );
    %  refined quasistatic elements
    [ y1, y2 ] = interp2( obj, green, data.pos1, data.pos2, 'ind', [ i, q, k ] );
    
    %  unit vectors in radial, azimuthal and z-direction
    [ er, et, ez ] = deal( data.er, data.et, data.ez );
    %  shape functions
    [ f, fp ] = deal( data.f, data.fp );
    fr = dot( er, f, k );
    ft = dot( et, f, k );
    fz = dot( ez, f, k );
    
    %  multiply SL terms with shape functions
    yy = struct;
    yy.pp  = y2.zz      * fp;
    yy.zp1 = y1.z1 * ez * fp;
    yy.zp2 = y1.z2 * ez * fp;
    yy.pz1 = y2.z1      * fz;
    yy.pz2 = y2.z2      * fz;
    yy.zz  = y1.zz * ez * fz;
    yy.ss  = y1.s * ( er * fr + et * ft );
    %  multiply DL terms with shape functions
    yy.pr = y1.rz1 * et * fr ./ data.r - y2.rz1 * ft;
    yy.rp = y1.rz2 * et * fp;
    yy.zr = y1.r   * ez * ft;
    yy.rz = y1.r   * et * fz;
    
    %  perform integration
    for name = convertCharsToStrings( fieldnames( yy ) ) .'
      yy.( name ) = double( sum( yy.( name ), q ), [ i, a, k ] );
    end
    %  indices
    yy.i1 = obj.yout( it ).i1( pt.i1 );
    yy.i2 = obj.yout( it ).i2( pt.i1 );
    %  add to output
    yout{ end + 1 } = yy;
  end
end
    
%  vectorize output
yout = horzcat( yout{ : } );
%  concatenate refined elements
obj.yout = struct;
for name = convertCharsToStrings( fieldnames( yout ) ) .'
  obj.yout.( name ) = cat( 1, yout.( name ) );
end
