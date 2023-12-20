function in = inner( obj, k0, varargin )
%  INNER - Inner product of Galerkin shape functions.
%
%  Usage for obj = feibelman.bemsolver :
%    in = inner( tau, k0, inout, PropertyPairs )
%  Input
%    tau      :  boundary elements
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    rules  :  quadrature rules
%  Output
%    in.I1    :  inner product <f,f>
%    in.I2    :  inner product <f,cross(nvec,f)>
%    in.d1    :  inner product <div(f),d1*div(f)>
%    in.d2    :  inner product <f,d2*cross(nvec,f)>

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'rules', quadboundary.rules );
%  parse input
parse( p, varargin{ : } );

%  allocate output
tau = obj.tau;
n = ndof( tau );
[ in.I1, in.I2, in.d11, in.d12, in.d21, in.d22 ] = deal( sparse( n, n ) );

for pt = quadboundary( tau, p.Results.rules )
  %  quadrature positions, weights, and shape elements
  [ ~, w, f, fp ] = eval( pt );
  %  dummy indices for tensor class
  [ i, q, a1, a2, k ] = deal( 1, 2, 3, 4, 5 );
  w = tensor( w, [ i, q ] );
  %  normal vector
  nvec = tensor( vertcat( pt.tau.nvec ), [ i, k ] );
    
  %  permittivity function
  mat = pt.mat( pt.inout( 1 ) );  eps1 = mat.eps( k0 );
  mat = pt.mat( pt.inout( 2 ) );  eps2 = mat.eps( k0 );
  %  Feibelman parameters
  [ d1, d2 ] = eval( obj.param, pt.tau, k0 );
  %  convert to tensor class
  d1 = tensor( d1, i );
  d2 = tensor( d2, i );
    
  %  integrands
  I1 = dot( tensor( f, [ i, q, a1, k ] ), tensor( f, [ i, q, a2, k ] ), k );
  I2 = dot( tensor( f, [ i, q, a1, k ] ),  ...
                             cross( nvec, tensor( f, [ i, q, a2, k ] ), k ), k );
  I3 = tensor( fp, [ i, q, a1 ] ) * tensor( fp, [ i, q, a2 ] );      
  %  perform integration
  d11 = 1i / ( k0 * eps1 ) * double( sum( w * d1 * I3, q ), [ i, a1, a2 ] );
  d12 = 1i / ( k0 * eps2 ) * double( sum( w * d1 * I3, q ), [ i, a1, a2 ] );
  d21 = 1i * ( k0 * eps1 ) * double( sum( w * d2 * I2, q ), [ i, a1, a2 ] );
  d22 = 1i * ( k0 * eps2 ) * double( sum( w * d2 * I2, q ), [ i, a1, a2 ] );
  
  I1 = double( sum( w * I1, q ), [ i, a1, a2 ] );
  I2 = double( sum( w * I2, q ), [ i, a1, a2 ] );
     
  %  convert from boundary elements to global degrees of freedom
  for it = 1 : numel( pt.tau )
    for a1 = 1 : 3
    for a2 = 1 : 3
      [ nu1, nu2 ] = deal( pt.tau( it ).nu( a1 ), pt.tau( it ).nu( a2 ) );
      in.I1( nu1, nu2 ) = in.I1( nu1, nu2 ) + I1( it, a1, a2 );
      in.I2( nu1, nu2 ) = in.I2( nu1, nu2 ) + I2( it, a1, a2 );
      in.d11( nu1, nu2 ) = in.d11( nu1, nu2 ) + d11( it, a1, a2 );
      in.d12( nu1, nu2 ) = in.d12( nu1, nu2 ) + d12( it, a1, a2 );
      in.d21( nu1, nu2 ) = in.d21( nu1, nu2 ) + d21( it, a1, a2 );
      in.d22( nu1, nu2 ) = in.d22( nu1, nu2 ) + d22( it, a1, a2 );      
    end
    end
  end
end
