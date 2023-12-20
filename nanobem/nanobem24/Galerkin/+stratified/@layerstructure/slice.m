function pts = slice( obj, varargin )
%  SLICE - Group positions into unique layer slices.
%
%  Usage for obj = stratified.layerstructure :
%    pts = slice( obj, pos1 )
%    pts = slice( obj, pos1, pos2 )
%  Input
%    pos1   :  first  set of positions
%    pos2   :  second set of positions
%  Output
%    pts    :  grouped positions or position pairs

switch numel( varargin )
  case 1
    %  layer indices
    pos = varargin{ 1 };
    ind = indlayer( obj, pos );
    %  allocate output
    i1 = unique( ind );
    pts = cell( 1, numel( i1 ) );
    %  loop over unique layer indices
    for it = 1 : numel( i1 )  
      k = i1( it );
      pts{ it } = struct(  ...
        'i1', i1( it ), 'pos', pos( ind == k, : ), 'ind', find( it == k ) );
    end
    
  otherwise
    %  layer indices 
    [ pos1, pos2 ] = deal( varargin{ : } );
    ind1 = indlayer( obj, pos1 );  i1 = unique( ind1 );
    ind2 = indlayer( obj, pos2 );  i2 = unique( ind2 );

    %  allocate output
    pts = cell( obj.n + 1 );
    %  loop over unique layer indces
    for it1 = 1 : numel( i1 );  k1 = i1( it1 );
    for it2 = 1 : numel( i2 );  k2 = i2( it2 );
      pts{ k1, k2 }.i1 = k1;
      pts{ k1, k2 }.i2 = k2;
      pts{ k1, k2 }.pos1 = pos1( ind1 == k1, : );
      pts{ k1, k2 }.pos2 = pos2( ind2 == k2, : );
      %  position indices
      pts{ k1, k2 }.ind1 = find( ind1 == k1 );
      pts{ k1, k2 }.ind2 = find( ind2 == k2 );
    end
    end    
end

%  make output vector
pts = horzcat( pts{ : } );
