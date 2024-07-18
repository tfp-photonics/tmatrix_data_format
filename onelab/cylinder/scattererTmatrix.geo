///////////////////////////////
// Author : Guillaume Demesy //
// scattererTmatrix.geo      //
// modified by R. Vallée for real-valued diagonally anisotropic material and lossy material 
///////////////////////////////

Include "scattererTmatrix_data.geo";
SetFactory("OpenCASCADE");

In_n  =  Sqrt[Fabs[epsr_In_re]];

paramaille_pml = paramaille/1.1;

If (flag_lossy==1)
  Out_lc        = 3*rbb/paramaille;
  PML_lc        = 3*rbb/paramaille_pml;
  In_lc         = 3*rbb/(paramaille*refine_scat);
Else
  Out_lc        = lambda_bg/paramaille;
  PML_lc        = lambda_bg/paramaille_pml;
  In_lc         = lambda0/(paramaille*In_n*refine_scat);
EndIf

If(flag_shape==ELL)
  Sphere (1) = { 0,0,0,ell_rx};
  Dilate { { 0,0,0 }, { 1, ell_ry/ell_rx, ell_rz/ell_rx } } { Volume{1}; }
EndIf
If(flag_shape==PARALL)
  Box (1) = {-par_ax/2,-par_ay/2,-par_az/2,par_ax,par_ay,par_az}; 
EndIf
If(flag_shape==CYL)
  Cylinder (1) = {0,0,-cyl_h/2,0,0,cyl_h,cyl_rx};
  Dilate { { 0,0,0 }, { 1 , cyl_ry/cyl_rx , 1} } { Volume{1}; }
EndIf
If (flag_shape==CONE)
  Cone (1) = {0,0,-cone_h/2,0,0,cone_h,cone_rx,0};
  Dilate { { 0,0,0 }, { 1 , cone_ry/cone_rx , 1} } { Volume{1}; }
EndIf
If (flag_shape==TOR)
  U (1) = {0,0,0,tor_r1,tor_r2x,tor_angle*Pi/180};
  Dilate { { 0,0,0 }, { 1 , 1 , tor_r2z/tor_r2x} } { Volume{1}; }
EndIf

Sphere (2) = { 0,0,0,r_pml_in}; 
Sphere (3) = { 0,0,0,r_pml_out}; 

Coherence;

Physical Volume("Scatterer" ,1) = {1};
Physical Volume("Background",2) = {2};
Physical Volume("PML"       ,3) = {3};

If(flag_shape==ELL)
  Physical Surface("SurfInt",10) = {2};
  Physical Surface("SurfPML",20) = {3};
EndIf
If(flag_shape==PARALL)
  Physical Surface("SurfInt",10) = {7};
  Physical Surface("SurfPML",20) = {8};
EndIf
If(flag_shape==CYL)
  Physical Surface("SurfInt",10) = {4};
  Physical Surface("SurfPML",20) = {5};
EndIf
If(flag_shape==CONE)
  Physical Surface("SurfInt",10) = {3};
  Physical Surface("SurfPML",20) = {4};
EndIf
If(flag_shape==TOR)
  Physical Surface("SurfInt",10) = {4};
  Physical Surface("SurfPML",20) = {5};
EndIf

Characteristic Length{PointsOf{Physical Volume{3};}} = PML_lc;
Characteristic Length{PointsOf{Physical Volume{2};}} = Out_lc;
Characteristic Length{PointsOf{Physical Volume{1};}} = In_lc;

Solver.AutoMesh=2;
Geometry.Points = 1;
Mesh.VolumeEdges = 0;

Mesh.ElementOrder = 1;
// Mesh.ElementOrder = 2;
// Mesh.HighOrderOptimize = 1;
