function y = compress( ~, y )
%  COMPRESS - Compress structure to single ODE vector.

y = cellfun( @( name ) reshape( y.( name ), [], 1 ), fieldnames( y ), 'uniform', 0 );
y = vertcat( y{ : } );
