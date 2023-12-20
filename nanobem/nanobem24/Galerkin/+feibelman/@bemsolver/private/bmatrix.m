function [ b1, b2, I ] = bmatrix( obj, k0 )
%  BMATRIX - Boundary matrices for BEM solver.
%
%  Usage for obj = feibelman.bemsolver :
%    [ b1, b2, I ] = bmatrix( obj, k0 )
%  Input
%    k0     :  wavenumber of light in vacuum
%  Output
%    b1     :  boundary matrix at particle inside
%    b2     :  boundary matrix at particle outside
%    I      :  identity matrix for Calderon indentity

%  number of global degrees of freedom
tau = obj.tau;
n = ndof( tau );
%  inner products of Galerkin shape functions
in = inner( obj, k0 );

%  identity matrix
nul = zeros( n );
I = - [ in.I2, nul; nul, in.I2 ];

%  boundary condition matrices
I1 = in.I1;
b1 = [ I1, ( in.I2 / I1 ) * in.d11; in.d21, I1 ];
b2 = [ I1, ( in.I2 / I1 ) * in.d12; in.d22, I1 ]; 
