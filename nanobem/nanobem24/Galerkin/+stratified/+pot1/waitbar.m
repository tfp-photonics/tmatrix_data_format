function y = waitbar( mode, varargin )
%  WAITBAR - Show waitbar during potential initialization.

switch mode
  case 'init'
    %  extract input
    [ x, y.waitbar, y.name ] = deal( varargin{ : } );
    %  number of integration points
    ind = arrayfun( @( x ) ~isempty( x.yout ), x, 'uniform', 1 );
    n = arrayfun(  ...
      @( x ) stratified.pot1.npts( x.yout ), x( ind ), 'uniform', 0 );
    [ y.n, y.it ] = deal( horzcat( n{ : } ), 0 );
    %  show waitbar during initialization ?
    if y.waitbar,  multiWaitbar( y.name, 'Color', 'r' );  end

  case 'close'
    %  extract input and close waitbar
    y = varargin{ 1 };
    if y.waitbar,  multiWaitbar( y.name, 'close' );  end
    
  otherwise
    %  extract input and update waitbar
    y = varargin{ 1 };
    y.it = y.it + 1;
    if y.waitbar
      multiWaitbar( y.name, sum( y.n( 1 : y.it ) ) / sum( y.n ) );  
    end
end