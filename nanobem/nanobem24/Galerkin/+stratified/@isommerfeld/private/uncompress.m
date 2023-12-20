function yout = uncompress( ~, y, yout )
%  UNCOMPRESS - Uncompress structure from single ODE vector.
%
%  Usage for obj = stratified.isommerfeld :
%    yout = uncompress( obj, y, yout )
%  Input
%    y      :  ODE vector
%    yout   :  template for output structure

%  size of output arrays
names = fieldnames( yout );
siz = cellfun( @( name ) size( yout.( name ) ), names, 'uniform', 0 );

%  split ODE vector and set output structure
y = mat2cell( y( : ), cellfun( @( siz ) prod( siz ), siz, 'uniform', 1 ) );
for it = 1 : numel( y )
  yout.( names{ it } ) = reshape( y{ it }, siz{ it } );
end
