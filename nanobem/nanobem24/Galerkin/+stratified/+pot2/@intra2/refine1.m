function obj = refine1( obj )
%  REFINE1 - Refine quasistatic terms.
%
%  Usage for obj = stratified.pot2.intra2 :
%    obj = refine1( obj )

%  Green function object and allocate output
green = stratified.tab.intra2( obj.layer, obj.i1, [], [], [] );
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
    
    yy = [ struct, struct ];
    %  loop over upper and lower interface
    for m = 1 : 2
      %  multiply SL terms with shape functions
      yy( m ).x  = y1( m ).z  * ez * fp - y2( m ).z * fz;
      yy( m ).pp = y2( m ).zz      * fp;
      yy( m ).zz = y1( m ).zz * ez * fz;
      yy( m ).ss = y1( m ).s * ( er * fr + et * ft );
      %  multiply DL terms with shape functions
      yy( m ).pr = y1( m ).rz * et * fr ./ data.r - y2( m ).rz * ft;
      yy( m ).rp = y1( m ).rz * et * fp;
      yy( m ).zr = y1( m ).r  * ez * ft;
      yy( m ).rz = y1( m ).r  * et * fz;
    
      %  perform integration
      for name = convertCharsToStrings( fieldnames( yy( m ) ) ) .'
        yy( m ).( name ) = double( sum( yy( m ).( name ), q ), [ i, a, k ] );
      end
    end
    %  indices
    [ yy( 1 ).i1, yy( 2 ).i1 ] = deal( obj.yout( it ).i1( pt.i1 ) );
    [ yy( 1 ).i2, yy( 2 ).i2 ] = deal( obj.yout( it ).i2( pt.i1 ) );    
    %  add to output
    yout{ end + 1, 1 } = yy( 1 );
    yout{ end + 1, 2 } = yy( 2 );
  end
end
    
%  vectorize output
yout1 = horzcat( yout{ :, 1 } );
yout2 = horzcat( yout{ :, 2 } );
%  concatenate refined elements
obj.yout = [ struct, struct ];
for name = convertCharsToStrings( fieldnames( yout1 ) ) .'
  obj.yout( 1 ).( name ) = cat( 1, yout1.( name ) );
  obj.yout( 2 ).( name ) = cat( 1, yout2.( name ) );
end
