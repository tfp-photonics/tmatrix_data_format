function [ e, h ] = planewave( obj, pol, dir, k0, pos, varargin )
%  PLANEWAVE - Propagate plane wave through layer structure.
%
%  Usage for obj = stratified.layerstructure :
%    [ e, h ] = planewave( obj, pol, dir, k0, pos, PropertyPairs )
%  Input
%    pol      :  light polarizations
%    dir      :  light propagation directions
%    k0       :  wavenumber of light in vacuum
%    pos      :  positions where fields are evaluated
%  PropertyName
%    primary  :  add primary fields ?
%  Output
%    e        :  electric field
%    h        :  magnetic field

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'primary', 1 );
%  parse input
parse( p, varargin{ : } );

%  number of light propagation directions
n = size( dir, 1 );
%  wavenumbers and impedances
k = arrayfun( @( x ) x.k( k0 ), obj.mat, 'uniform', 1 );
Z = arrayfun( @( x ) x.Z( k0 ), obj.mat, 'uniform', 1 );
%  material index of incoming light and first interface
switch sign( dir( 1, 3 ) )
  case 1
    [ i2, z2 ] = deal( 1, obj.z( 1 ) );
    Z1 = Z( 1 );
  otherwise
    [ i2, z2 ] = deal( obj.n + 1, obj.z( end ) );    
    Z1 = Z( end );
end
%  wavenumber of incoming light
k2 = k( i2 ) * dir;
kpar = sqrt( k2( :, 1 ) .^ 2 + k2( :, 2 ) .^ 2 );

%  normalization function 
norm = @( x ) bsxfun( @rdivide, x, sqrt( sum( abs( x ) .^ 2, 2 ) ) );
i1 = abs( dir( :, 3 ) ) == 1;
%  unit vectors for TE and TM decomposition
te = norm( cross( dir, repmat( [ 0, 0, 1 ], n, 1 ) ) );
te( i1, : ) = norm( pol( i1, : ) );
tm = cross( conj( te ), dir, 2 ); 
%  amplitudes for TE and TM waves
e2.te = dot( te, pol, 2 ) .* exp( 1i * k2( :, 3 ) * z2 );
e2.tm = dot( tm, pol, 2 ) .* exp( 1i * k2( :, 3 ) * z2 ) / Z1;

%  layer indices for observation points
i1 = indlayer( obj, pos );
%  allocate output
[ e, h ] = deal( zeros( [ size( pos ), n ] ) );

%  loop over unique indices
for it = unique( i1 .' )
  %  coefficients for secondary waves
  wave = secondary( obj, k0, kpar, it, i2 ); 
  ind = i1 == it; 
  
  %  loop over TE and TM polarizations
  for name = [ "te", "tm" ]
    a = horzcat( wave.( name ) ) .';
    kz = vertcat( wave.kz );
    kz = kz( :, it );
    %  loop over wave directions
    for d1 = 1 : size( a, 2 )
      %  propagation direction and distance to interface
      if d1 == 1 && it ~= 1
        dir1 = 1;
        z = pos( ind, 3 ) - obj.z( it - 1 );
      else
        dir1 = - 1;
        z = obj.z( it ) - pos( ind, 3 ); 
      end
      %  secondary fields
      xpar = k2( :, 1 : 2 ) * pos( ind, 1 : 2 ) .' + kz * z .';
      e1 = bsxfun( @times, exp( 1i * xpar ), a( :, d1 ) .* e2.( name ) );
      %  unit propagation vector
      dir1 = [ k2( :, 1 : 2 ), dir1 * kz ] / obj.mat( it ).k( k0 );  
      %  convert to tensor
      [ i, nu, k ] = deal( 1, 2, 3 );
      dir1 = tensor( dir1, [ nu, k ] );
      
      switch name
        case 'te'
          e1 = tensor( e1, [ nu, i ] ) * tensor( te, [ nu, k ] );          
          h1 = cross( dir1, e1, k ) / Z( it );
        case 'tm'
          h1 = tensor( e1, [ nu, i ] ) *  ...
               tensor( cross( dir, tm, 2 ), [ nu, k ] );             
          e1 = cross( h1, dir1, k ) * Z( it );
      end
      %  add up secondary fields
      e( ind, :, : ) = e( ind, :, : ) + double( e1, [ i, k, nu ] ); 
      h( ind, :, : ) = h( ind, :, : ) + double( h1, [ i, k, nu ] ); 
    end
  end
  %  add primary field ?
  if p.Results.primary && it == i2
    a = exp( 1i * pos( ind, : ) * k2 .' );
    e1 = tensor( a, [ i, nu ] ) * tensor( pol, [ nu, k ]  );
    h1 = tensor( a, [ i, nu ] ) *  ...
        tensor( cross( dir, pol, 2 ), [ nu, k ] ) / Z( it );
    e( ind, :, : ) = e( ind, :, : ) + double( e1, [ i, k, nu ] );
    h( ind, :, : ) = h( ind, :, : ) + double( h1, [ i, k, nu ] );
  end
end

