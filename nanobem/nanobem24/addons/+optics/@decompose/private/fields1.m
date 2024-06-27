function [ e, h ] = fields1( obj, pts, imat )
%  FIELDS1 - Fields for planewave decomposition and unbounded medium.
%
%  Usage for obj = optics.decompose :
%    [ e, h ] = fields1( obj, pts, imat )
%  Input
%    pts    :  points where electromagnetic fields are requested
%    imat   :  material index for unbounded medium
%  Output
%    e,h    :  electromagnetic fields

%  index to points with material index IMAT
ind = vertcat( pts.imat ) == imat;
%  allocate output
[ e, h ] = deal( zeros( numel( pts ), 3 ) );

%  evaluate fields
if nnz( ind )
  %  wavelength and impedance in medium
  mat = pts( find( ind, 1 ) ).mat( imat );
  [ k1, Z1 ] = deal( mat.k( obj.k0 ), mat.Z( obj.k0 ) );
  
  %  electromagnetic fields
  fac = exp( 1i * k1 * vertcat( pts( ind ).pos ) * obj.dir .' );
  e( ind, : ) = fac * obj.efield;
  h( ind, : ) = fac * cross( obj.dir, obj.efield, 2 ) / Z1;
end
