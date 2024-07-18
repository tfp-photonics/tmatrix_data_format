///////////////////////////////
// Author : Guillaume Demesy //
// scattererTmatrix_data.geo //
// modified by R. Vallée for real-valued diagonally anisotropic material and lossy material 
///////////////////////////////

nm       = 1.;
epsilon0 = 8.854187817e-3*nm;
mu0      = 400.*Pi*nm;
cel      = 1.0/(Sqrt[epsilon0 * mu0]);
deg2rad  = Pi/180.;

pp0        = "1Geometry/0";
pp1        = "0Study Type/0";
pp2        = "3Electromagnetic parameters/0";
pp3        = "4Mesh size and PMLs parameters/0";
pp4        = "5Postpro options/0";
pp5        = "6Post plot options/0";
close_menu = 0;
colorro    = "LightGrey";
colorppOK   = "Ivory";
colorppWA   = "LightSalmon1";
colorppNO   = "LightSalmon4";

ELL    = 0;
PARALL = 1;
CYL    = 2;
CONE   = 3;
TOR    = 4;

RES_PW    = 0;
RES_TMAT  = 1;
RES_GREEN = 2;
RES_QNM   = 3;

DefineConstant[
  flag_shape = {PARALL , Name StrCat[pp0, "0Scatterer shape"], 
    Choices {ELL="ellispoid",PARALL="parallelepiped",CYL="cylinder",CONE="cone",TOR="split U"},Closed 0}];

If (flag_shape==ELL)
  DefineConstant[
    ell_rx = {150.0 , Name StrCat[pp0  , "1ellipsoid X-radius [nm]"], Highlight Str[colorppOK] , Closed 0},
    ell_ry = {150.0 , Name StrCat[pp0  , "2ellipsoid Y-radius [nm]"], Highlight Str[colorppOK] , Closed 0},
    ell_rz = {150.0 , Name StrCat[pp0  , "3ellipsoid Z-radius [nm]"], Highlight Str[colorppOK] , Closed 0}
  ];
  // FIXME gmsh Max?
  If (ell_rx>ell_ry)
    rbb_tmp=ell_rx;
  Else rbb_tmp=ell_ry;
  EndIf
  If (ell_rz>rbb_tmp)
    rbb=ell_rz;
  Else rbb=rbb_tmp;
  EndIf
EndIf
If (flag_shape==PARALL)
  DefineConstant[
    par_ax = {250 , Name StrCat[pp0  , "1cube X-edge size [nm]"], Highlight Str[colorppOK]  , Closed 0},
    par_ay = {200 , Name StrCat[pp0  , "2cube Y-edge size [nm]"], Highlight Str[colorppOK]  , Closed 0},
    par_az = {150 , Name StrCat[pp0  , "3cube Z-edge size [nm]"], Highlight Str[colorppOK]  , Closed 0}];
  rbb = 0.5*Sqrt[par_ax^2+par_ay^2+par_az^2];
EndIf
If (flag_shape==CYL)
  DefineConstant[
    cyl_rx = {250 , Name StrCat[pp0  , "1cylinder X-radius [nm]"], Highlight Str[colorppOK]  , Closed 0},
    cyl_ry = {250 , Name StrCat[pp0  , "2cylinder Y-radius [nm]"], Highlight Str[colorppOK]  , Closed 0},
    cyl_h  = {300 , Name StrCat[pp0  , "3cylinder height [nm]"]  , Highlight Str[colorppOK]  , Closed 0}];
  // FIXME gmsh Max?
  If (cyl_rx>cyl_ry)
    rbb_tmp=cyl_rx;
  Else rbb_tmp=cyl_ry;
  EndIf
  rbb = 0.5*Sqrt[rbb_tmp^2+cyl_h^2];
EndIf
If (flag_shape==CONE)
  DefineConstant[
    cone_rx = {300 , Name StrCat[pp0  , "1cone basis X-radius [nm]"], Highlight Str[colorppOK]  , Closed 0},
    cone_ry = {300 , Name StrCat[pp0  , "1cone basis Y-radius [nm]"], Highlight Str[colorppOK]  , Closed 0},
    cone_h  = {300 , Name StrCat[pp0  , "2cone height [nm]"]        , Highlight Str[colorppOK]  , Closed 0}];
  // FIXME gmsh Max?
  If (cone_rx>cone_ry)
    rbb_tmp=cone_rx;
  Else rbb_tmp=cone_ry;
  EndIf    
  If (2.0*cone_h/3.0>Sqrt[rbb_tmp^2+cone_h^2/9.0])
    rbb=2.0*cone_h/3.0;
  Else rbb=Sqrt[rbb_tmp^2+cone_h^2/9.0];
  EndIf
EndIf
If (flag_shape==TOR)
  DefineConstant[
    tor_r1    = { 300 , Name StrCat[pp0 , "1U radius 1  [nm]"], Highlight Str[colorppOK] , Closed 0},
    tor_r2x   = { 100 , Name StrCat[pp0 , "2U radius 2x [nm]"], Highlight Str[colorppOK] , Closed 0},
    tor_r2z   = { 50  , Name StrCat[pp0 , "3U radius 2z [nm]"], Highlight Str[colorppOK] , Closed 0},
    tor_angle = { 340 , Name StrCat[pp0 , "4U angle [deg]"]   , Highlight Str[colorppOK] , Closed 0, Min 5, Max 355}];
  rbb = tor_r1+tor_r2x;
EndIf
DefineConstant[
  rot_theta = {0 , Name StrCat[pp0  , "5rotate scatterer (polar) [deg]"] , Highlight Str[colorppOK]  , Closed 0, Min 0, Max 180},
  rot_phi   = {0 , Name StrCat[pp0  , "6rotate scatterer (azimut) [deg]"], Highlight Str[colorppOK]  , Closed 0, Min 0, Max 360}];

flag_study = RES_TMAT;

DefineConstant[
  epsr_In_re  = { 9, Name StrCat[pp2  , "0scatterer permittivity (real) []"] , Highlight Str[colorppOK] , Closed 0},
  epsr_In_im  = { 0., Name StrCat[pp2  , "1scatterer permittivity (imag) []"] , Highlight Str[colorppOK] , Closed 0},
  epsr_Out_re = { 1., Name StrCat[pp2  , "2background permittivity (real) []"], Highlight Str[colorppOK] , Closed 0},
  epsr_Out_im = { 0., Name StrCat[pp2  , "3background permittivity (imag) []"], Highlight Str[colorppOK] , Closed 0}
  ];

DefineConstant[
  flag_anisotropy = {0, Choices{0,1}, Name StrCat[pp2, "anisotropic structure?"]}
];

//only take into account real-valued diagonal anisotropy; is obvious to extend
If (flag_anisotropy==1)
  DefineConstant[
    epsr_In_re_xx  = { 3.24, Name StrCat[pp2  , "0xx scatterer permittivity (re)"]},
    epsr_In_re_yy  = { 4.0, Name StrCat[pp2  , "0yy scatterer permittivity (re)"]},
    epsr_In_re_zz  = { 4.84, Name StrCat[pp2  , "0zz scatterer permittivity (re)"]}
  ];
EndIf

DefineConstant[
  lambda0 = {600 , Name StrCat[pp2  , "4wavelength [nm]"] , Highlight Str[colorppOK] , Closed 0, GmshOption "Reset", Autocheck 0}
];

lambda_bg = lambda0/Sqrt[epsr_Out_re];
  
DefineConstant[
  n_max = {3 , Name StrCat[pp2  , "8n_max integer"], Highlight Str[colorppOK]  , Closed 0, Min 0, Max 5},
  siwt  = {0 , Name StrCat[pp2  , "9Time sign e^(+|-iwt)"], Choices {0="e^(-iwt)",2="e^(+iwt)"}  , Closed 0}
];

flag_cartpml = 0;

DefineConstant[
  flag_lossy = {0, Choices{0,1}, Name StrCat[pp2, "lossy structure?"]}
];

If (flag_lossy==1)
  DefineConstant[
    space2pml    = {1.5*rbb , Name StrCat[pp3  , "0space around scatterer [nm]"] , Min (rbb*0.2), Max (rbb*5)},
    pml_size     = {3*rbb , Name StrCat[pp3  , "1PML thickness [nm]"] , Highlight Str[colorppWA]  , Closed 0, Min (lambda_bg/10), Max (3*lambda_bg)}    
  ];
Else
  DefineConstant[
    space2pml    = {lambda_bg/2 , Name StrCat[pp3  , "0space around scatterer [nm]"] , Highlight Str[colorppNO] , Closed 1,Min (lambda_bg/10.), Max (3*lambda_bg)},
    pml_size     = {lambda_bg , Name StrCat[pp3  , "1PML thickness [nm]"] , Highlight Str[colorppWA]  , Closed 0, Min (lambda_bg/10), Max (3*lambda_bg)}
  ];
EndIf

DefineConstant[
  refine_scat  = {1. , Name StrCat[pp3  , "3scatterer mesh refinement"], Highlight Str[colorppWA] , Closed 0}, 
  paramaille   = {5. , Name StrCat[pp3  , "2mesh size"], Highlight Str[colorppWA]  , Closed 0},
  is_FEM_o2    = {1 , Name StrCat[pp3   , "4Interpolation order "] , Choices {0="order 1",1="order 2"}, Closed 0}
];

DefineConstant[
  flag_plotcuts = {0, Choices{0,1}, Name StrCat[pp5, "Plot radial cuts?"]}
];

// FIXME conditional for normalization
If (flag_shape==ELL)
  ell_rx = ell_rx*nm;
  ell_ry = ell_ry*nm;
  ell_rz = ell_rz*nm;
EndIf
If (flag_shape==PARALL)
    par_ax = par_ax*nm;
    par_ay = par_ay*nm;
    par_az = par_az*nm;
EndIf
If (flag_shape==CYL)
  cyl_rx = cyl_rx*nm;
  cyl_ry = cyl_ry*nm;
  cyl_h  = cyl_h *nm;
EndIf
If (flag_shape==CONE)
  cone_rx = cone_rx*nm;
  cone_ry = cone_ry*nm;
  cone_h = cone_h*nm;
EndIf
If (flag_shape==TOR)
  tor_r1 = tor_r1*nm;
  tor_r2x = tor_r2x*nm;
  tor_r2z = tor_r2z*nm;
EndIf

lambda0   = lambda0*nm;
lambda_bg = lambda_bg*nm;
pml_size  = pml_size*nm;
space2pml = space2pml*nm;
rbb       = rbb*nm;

r_pml_in  = space2pml+rbb;
r_pml_out = space2pml+rbb+pml_size;


p_max = n_max*n_max+2*n_max;
siwt=siwt-1;

sph_scan = 0.00001;

npts_plot_theta = 25;
npts_plot_phi   = 50;
