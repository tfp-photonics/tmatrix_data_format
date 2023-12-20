function f = fields( obj, k0, kpar, z1, z2, varargin )
%  EFIELD - Fields at positions Z1 for source at position Z2.
%
%  Usage for obj = stratified.layerstructure :
%    f = fields( obj, k0, kpar, z1, z2, PropertyPairs )
%  Input
%    k0     :  wavenumber in vacuum
%    kpar   :  parallel momentum
%    z1     :  observation points
%    z2     :  source point, must be scalar
%  PropertyName
%     primary   :  add primary fields ?
%  Output
%    f          :  electric (TE) or magnetic (TM) field

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'primary', 1 );
%  parse input
parse( p, varargin{ : } );

%  function mainly for testing, use single source point only
assert( isscalar( z2 ) );
%  layer indices
i1 = indlayer( obj, z1( : ) * [ 0, 0, 1 ] );
i2 = indlayer( obj, [ 0, 0, z2 ] );
%  allocate output
z1 = z1( : );
[ f.te, f.tm ] = deal( 0 * z1 );
%  wave direction
dir = [ 1, -1 ];

%  loop over unique indices
for it = unique( i1 .' )
  %  coefficients for secondary waves
  wave = secondary( obj, k0, kpar, it, i2 );
  %  wavenumbers 
  kz1 = wave.kz( it );
  kz2 = wave.kz( i2 );
  ind = i1 == it; 
  %  loop over TE and TM polarizations
  for name = [ "te", "tm" ]
    a = wave.( name );
    %  loop over wave directions
    for d1 = 1 : size( a, 1 )
    for d2 = 1 : size( a, 2 )
      dir1 = dir( d1 );  if it == 1,          dir1 = -1;  end
      dir2 = dir( d2 );  if i2 == obj.n + 1,  dir2 = -1;  end
      %  distance to interfaces
      Z1 = fun( obj, z1( ind ), it, dir1 );
      Z2 = fun( obj, z2, i2, - dir2 );
      %  add up secondary fields
      f.( name )( ind ) =  ...
      f.( name )( ind ) + a( d1, d2 ) * exp( 1i * ( kz1 * Z1 + kz2 * Z2 ) ); 
    end
    end
    %  add primary field ?
    if p.Results.primary && it == i2
      f.( name )( ind ) =  ...
      f.( name )( ind ) + exp( 1i * kz1 * abs( z1( ind ) - z2 ) );
    end
  end
end


function z = fun( obj, z, it, dir )
%  FUN - Distance from interfaces.
switch dir
  case 1
    z = z - obj.z( it - 1 ); 
  otherwise
    z = obj.z( it ) - z; 
end
