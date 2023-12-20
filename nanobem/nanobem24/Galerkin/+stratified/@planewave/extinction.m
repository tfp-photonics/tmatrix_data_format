function cext = extinction( obj, sol, key )
%  EXTINCTION - Extinction cross section for planewave excitation.
%    See Lytle et al., Phys. Rev. E 71, 056610 (2005).
%
%  Usage for obj = stratified.planewave :
%    cext = extinction( obj, sol, key )
%  Input
%    sol    :  solution of BEM equations
%    key   :  'refl', 'trans' or 'both'
%  Output
%    cext   :  extinction cross section

if ~exist( 'key', 'var' ),  key = 'both';  end

%  reflected and transmitted incoming far-fields
[ er, ~, kr ] = farfields( obj, sol.k0, 'refl'  );
[ et, ~, kt ] = farfields( obj, sol.k0, 'trans' );
%  allocate output for extinction cross section
cext = zeros( 1, size( obj.dir, 1 ) );

%  loop over propagation directions
for it = 1 : size( obj.dir, 1 )
  
  %  scattered farfield for reflection
  esr = farfields( sol, kr( it, : ) / norm( kr( it, : ) ) );
  %  scattered farfields for transmission
  if abs( imag( kt( it, 3 ) ) ) > 1e-10 || norm( kt ) < 1e-10
    %  evanescent field or absorbing background medium
    est = [ 0, 0, 0 ]; 
  else
    est = farfields( sol, kt( it, : ) / norm( kt( it, : ) ) );  
  end
  
  %  extinction of reflected and transmitted beam
  extr = 4 * pi / norm( kr( it, : ) ) * imag( dot( er( it, : ), esr ) );
  extt = 4 * pi / norm( kr( it, : ) ) * imag( dot( et( it, : ), est ) );
  %  set output
  switch key
    case 'refl'
      cext( it ) = extr;
    case 'trans'
      cext( it ) = extt;
    otherwise
      cext( it ) = extr + extt;
  end
end
