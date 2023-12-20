function n = npts( x )
%  NPTS - Number of quadrature points.

switch class( x )
  case 'quadduffy'
    n = x.npts;
  case 'quadboundary'
    n = x( 1 ).npts * numel( x( 2 ).quad.w );
  case 'struct'
    n = arrayfun( @( x ) stratified.pot1.npts( x.pts ), x, 'uniform', 1 );
end
