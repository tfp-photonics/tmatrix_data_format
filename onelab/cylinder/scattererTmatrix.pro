///////////////////////////////
// Author : Guillaume Demesy //
// scattererTmatrix.pro      //
// modified by R. Vallée for real-valued diagonally anisotropic material and lossy material 
///////////////////////////////

Include "scattererTmatrix_data.geo";
myDir = StrCat["resT",StrCat[Sprintf["_wl_%g",lambda0],"/"]];

Group {
  Scat_In        = Region[1];
  Scat_Out       = Region[2];
  Domain         = Region[{Scat_In,Scat_Out}];
  PMLs           = Region[3];
  All_domains    = Region[{Scat_In,Scat_Out,PMLs}];
  SurfInt        = Region[{10}];
  SurfDirichlet  = Region[{20}];
}

Function{
  I[]         = Complex[0,1];
  avoid_sing  = 1.e-14;
  mu0         = 4*Pi*100.0*nm;
  epsilon0    = 8.854187817e-3*nm;
  cel         = 1.0/(Sqrt[epsilon0 * mu0]);
  k0          = 2.0*Pi/lambda0;
  k_Out       = 2.0*Pi*Sqrt[epsr_Out_re]/lambda0;
  R3D[]   = Sqrt[X[]^2+Y[]^2+Z[]^2];
  R2D[]   = Sqrt[X[]^2+Y[]^2];
  cos_theta[] = Z[]/(R3D[]+avoid_sing);
  theta[]     = Acos[cos_theta[]];
  phi[]       = Atan2[-Y[],-X[]]+Pi;
  cart2sph[]  =   Tensor[X[]/R3D[], Y[]/R3D[], Z[]/R3D[],       
                         X[]*Z[]/(R3D[]*R2D[]), Y[]*Z[]/(R3D[]*R2D[]), -(X[]^2+Y[]^2)/(R3D[]*R2D[]),
                        -Y[]/R2D[], X[]/R2D[],0];
  R[] = Vector[X[]/R3D[], Y[]/R3D[], Z[]/R3D[]];

  a_pml =  1.;
  b_pml = -siwt;

  sr[]         = Complex[a_pml,b_pml];
  sphi[]       = Complex[1,0];
  stheta[]     = Complex[1,0];
  r_tild[]     = r_pml_in + (R3D[] - r_pml_in) * sr[];
  theta_tild[] = theta[];
  pml_tens_temp1[]  = TensorDiag[(r_tild[]/R3D[])^2 * sphi[]*stheta[]/sr[],
                                 (r_tild[]/R3D[])   * sr[]*stheta[]/sphi[],
                                 (r_tild[]/R3D[])   * sphi[]*sr[]/stheta[]];
  pml_tens_temp2[]  = Rotate[pml_tens_temp1[],0,-theta[]-Pi/2,0];
  pml_tens[]        = Rotate[pml_tens_temp2[],0,0,-phi[]];

  epsilonr_In[]       = Complex[epsr_In_re  , epsr_In_im];
  epsilonr_Out[]      = Complex[epsr_Out_re , epsr_Out_im];

  If (flag_anisotropy==1)
    epsilonr[Scat_In]   = TensorDiag[epsr_In_re_xx,epsr_In_re_yy,epsr_In_re_zz];
  Else
    epsilonr[Scat_In]   = epsilonr_In[]  * TensorDiag[1.,1.,1.];
  EndIf
  
  epsilonr[Scat_Out]  = epsilonr_Out[] * TensorDiag[1.,1.,1.];
  epsilonr[PMLs]      = epsilonr_Out[] * pml_tens[];
                          
  epsilonr1[Scat_In]  = epsilonr_Out[] * TensorDiag[1.,1.,1.];
  epsilonr1[Scat_Out] = epsilonr_Out[] * TensorDiag[1.,1.,1.];                
  epsilonr1[PMLs]     = epsilonr_Out[] * pml_tens[];
                          
  mur[Scat_In]        = TensorDiag[1.,1.,1.];
  mur[Scat_Out]       = TensorDiag[1.,1.,1.];
  mur[PMLs]           = pml_tens[];

  For pe In {1:p_max}
    ne = Floor[Sqrt[pe]];
    me = ne*(ne+1) - Floor[pe];
    Mnm_source~{pe}[] = Mnm[1,ne,me,XYZ[],k_Out];
    Nnm_source~{pe}[] = Nnm[1,ne,me,XYZ[],k_Out];
    Mnm_out~{pe}[]    = Mnm[3,ne,me,XYZ[],k_Out];
    Nnm_out~{pe}[]    = Nnm[3,ne,me,XYZ[],k_Out];
    Xnm~{pe}[]        = Xnm[ne,me,X[],Y[],Z[]];
    Ynm~{pe}[]        = Ynm[ne,me,X[],Y[],Z[]];
    Znm~{pe}[]        = Znm[ne,me,X[],Y[],Z[]];
    Ynm_r~{pe}[]      = Ynm[ne,me,X[],Y[],Z[]] * R[];
    source_M~{pe}[]   = k0^2*(epsilonr[]-epsilonr1[])*Mnm_source~{pe}[];
    source_N~{pe}[]   = k0^2*(epsilonr[]-epsilonr1[])*Nnm_source~{pe}[];
    SphHankelOutgoing_n~{pe}[]   = (JnSph[ne  ,k_Out*r_pml_in]- siwt * I[]*YnSph[ne  ,k_Out*r_pml_in]) ;
    SphHankelOutgoing_nm1~{pe}[] = (JnSph[ne-1,k_Out*r_pml_in]- siwt * I[]*YnSph[ne-1,k_Out*r_pml_in]) ;
    dRicattiBessel~{pe}[]        = (k_Out * r_pml_in * (SphHankelOutgoing_nm1~{pe}[]-((ne+1)/((k_Out*r_pml_in))) * SphHankelOutgoing_n~{pe}[]) + SphHankelOutgoing_n~{pe}[]);
    normalize_fhnm_X~{pe}[]      = 1/SphHankelOutgoing_n~{pe}[];
    normalize_fenm_Z~{pe}[]      = k_Out*r_pml_in/dRicattiBessel~{pe}[];
    normalize_fenm_Y~{pe}[]      = k_Out*r_pml_in/(SphHankelOutgoing_n~{pe}[]*Sqrt[ne*(ne+1)]);
  EndFor  
}

Constraint {
  {Name Dirichlet; Type Assign;
   Case {
     { Region SurfDirichlet; Value 0.; }
   }
  }
}

Jacobian {
  { Name JVol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name JSur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
}

Integration {
  { Name Int_1 ;
    Case { 
      { Type Gauss ;
        Case { 
          { GeoElement Point        ; NumberOfPoints   1 ; }
          { GeoElement Line2        ; NumberOfPoints   4 ; }
          { GeoElement Triangle     ; NumberOfPoints  12 ; }
          { GeoElement Triangle2    ; NumberOfPoints  12 ; }
          { GeoElement Tetrahedron  ; NumberOfPoints  17 ; }
          { GeoElement Tetrahedron2 ; NumberOfPoints  17 ; }
        }
      }
    }
  }
}

FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge;
        Support Region[{All_domains,SurfInt}]; Entity EdgesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Edge_2E;
        Support Region[{All_domains,SurfInt}]; Entity EdgesOf[All]; }
        If (is_FEM_o2==1)
          { Name sn3; NameOfCoef un3; Function BF_Edge_3F_b;
            Support Region[{All_domains,SurfInt}]; Entity FacetsOf[All]; }
          { Name sn4; NameOfCoef un4; Function BF_Edge_3F_c;
            Support Region[{All_domains,SurfInt}]; Entity FacetsOf[All]; }
          { Name sn5; NameOfCoef un5; Function BF_Edge_4E;
            Support Region[{All_domains,SurfInt}]; Entity EdgesOf[All]; }
         EndIf
    }
    Constraint {
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      If (is_FEM_o2==1)
        { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
        { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint Dirichlet; }  
        { NameOfCoef un5; EntityType EdgesOf  ; NameOfConstraint Dirichlet; }  
      EndIf
    }
  }
}

Formulation {
  {Name VPWMN_helmholtz_vector; Type FemEquation;
   Quantity {
     { Name u; Type Local; NameOfSpace Hcurl;}
    }
    Equation {
      Galerkin { [-1/mur[]*Dof{Curl u}             , {Curl u}]; In All_domains; Jacobian JVol; Integration Int_1;  }
      Galerkin { [k0^2*epsilonr[]*Dof{u} , {u}     ]; In All_domains; Jacobian JVol; Integration Int_1;  }
      Galerkin { [ k0^2*(epsilonr[]-epsilonr1[])*
                      ($isN ? Nnm[1,$NE,$ME,XYZ[],k_Out] : Mnm[1,$NE,$ME,XYZ[],k_Out]) 
                                                    , {u} ]; In Scat_In; Jacobian JVol; Integration Int_1;}      
    }
  }
}

Resolution {
  { Name res_VPWall_helmholtz_vector;
    System {
      { Name T; NameOfFormulation VPWMN_helmholtz_vector; Type ComplexValue; }
  }
  Operation {
    CreateDir[Str[myDir]];        
    For pe In {1:p_max}
      Evaluate[ $isN = 0 ];
      Evaluate[ $PE  = pe ];
      Evaluate[ $NE  = Floor[Sqrt[$PE]] ];
      Evaluate[ $ME  = $NE*($NE+1) - Floor[$PE] ];
      If (pe==1)
        Generate[T];
        Solve[T];
      EndIf
      GenerateRHS[T];
      SolveAgain[T];
      PostOperation[VPWM_postop~{pe}];
      Evaluate[$isN=1];
      GenerateRHS[T];
      SolveAgain[T];
      PostOperation[VPWN_postop~{pe}];
    EndFor
    }
  }
}

PostProcessing {
  For pe In {1:p_max}
    { Name VPWM_postpro~{pe}; NameOfFormulation VPWMN_helmholtz_vector; NameOfSystem T;
      Quantity {
        { Name E_scat           ; Value { Local { [{u}]; In All_domains; Jacobian JVol; } } }
        { Name H_scat           ; Value { Local { [siwt*I[]/(mur[]*mu0*k0*cel)*{Curl u}]; In All_domains; Jacobian JVol; } } }
        { Name Mnm_source~{pe}  ; Value { Local { [   Mnm_source~{pe}[] ]; In All_domains; Jacobian JVol; } } }

        For po In {1:p_max}
          { Name feM~{pe}~{po} ; Value { Integral { [ 1/(r_pml_in^2*epsilonr_Out[]) * normalize_fenm_Z~{po}[]*  {u}*Conj[Znm~{po}[]]  ]; In SurfInt; Integration Int_1 ; Jacobian JSur; } } }
          { Name fhM~{pe}~{po} ; Value { Integral { [ 1/(r_pml_in^2*epsilonr_Out[]) * normalize_fhnm_X~{po}[]*  {u}*Conj[Xnm~{po}[]]  ]; In SurfInt; Integration Int_1 ; Jacobian JSur; } } }
        EndFor
      }
    }
    { Name VPWN_postpro~{pe}; NameOfFormulation VPWMN_helmholtz_vector; NameOfSystem T;
      Quantity {
        { Name E_scat            ; Value { Local { [{u}]; In All_domains; Jacobian JVol; } } }
        { Name H_scat            ; Value { Local { [siwt*I[]/(mur[]*mu0*k0*cel)*{Curl u}]; In All_domains; Jacobian JVol; } } }
        { Name Nnm_source~{pe}   ; Value {    Local { [   Nnm_source~{pe}[] ]; In All_domains; Jacobian JVol; } } }
        For po In {1:p_max}
          { Name feN~{pe}~{po}  ; Value { Integral { [ 1/(r_pml_in^2*epsilonr_Out[]) *normalize_fenm_Z~{po}[]*  {u}*Conj[Znm~{po}[]]  ]; In SurfInt; Integration Int_1 ; Jacobian JSur; } } }
          { Name fhN~{pe}~{po}  ; Value { Integral { [ 1/(r_pml_in^2*epsilonr_Out[]) *normalize_fhnm_X~{po}[]*  {u}*Conj[Xnm~{po}[]]  ]; In SurfInt; Integration Int_1 ; Jacobian JSur; } } }
        EndFor
      }
    }
  EndFor
}


PostOperation {
  For pe In {1:p_max}
    {Name VPWM_postop~{pe}; NameOfPostProcessing VPWM_postpro~{pe} ;
      Operation {
        If (flag_plotcuts==1)
          Print [ E_scat , OnGrid
                {(r_pml_in-1*nm)*Sin[$B]*Cos[$C], (r_pml_in-1*nm)*Sin[$B]*Sin[$C], (r_pml_in-1*nm)*Cos[$B]} 
                {(r_pml_in-1*nm), {sph_scan : Pi-sph_scan+(Pi-2.0*sph_scan)/(10*(npts_plot_theta-1.0)) : (Pi-2.0*sph_scan)/(npts_plot_theta-1.0)}, 
                {sph_scan : 2.0*Pi-sph_scan+(2.0*Pi-2.0*sph_scan)/(10*(npts_plot_phi-1.0)) : (2.0*Pi-2.0*sph_scan)/(npts_plot_phi-1.0)} }, 
                File StrCat[myDir,StrCat["E_scat_onsphere_cart_M",Sprintf["%g.pos",pe]]],
                Name StrCat["E_scat_onsphere_cart_M",Sprintf["%g",pe]]];
        EndIf
        For po In {1:p_max}
          Print[feM~{pe}~{po}[SurfInt] , OnGlobal, Format Table, File StrCat[myDir,StrCat[StrCat["feM_",Sprintf["pe%g",pe]],Sprintf["po%g.dat",po]]]];
          Print[fhM~{pe}~{po}[SurfInt] , OnGlobal, Format Table, File StrCat[myDir,StrCat[StrCat["fhM_",Sprintf["pe%g",pe]],Sprintf["po%g.dat",po]]]];
        EndFor
      }
    }
    {Name VPWN_postop~{pe}; NameOfPostProcessing VPWN_postpro~{pe} ;
      Operation {
        If (flag_plotcuts==1)
          Print [ E_scat , OnGrid
                {(r_pml_in-1*nm)*Sin[$B]*Cos[$C], (r_pml_in-1*nm)*Sin[$B]*Sin[$C], (r_pml_in-1*nm)*Cos[$B]} 
                {(r_pml_in-1*nm), {sph_scan : Pi-sph_scan+(Pi-2.0*sph_scan)/(10*(npts_plot_theta-1.0)) : (Pi-2.0*sph_scan)/(npts_plot_theta-1.0)}, 
                {sph_scan : 2.0*Pi-sph_scan+(2.0*Pi-2.0*sph_scan)/(10*(npts_plot_phi-1.0)) : (2.0*Pi-2.0*sph_scan)/(npts_plot_phi-1.0)} }, 
                File StrCat[myDir,StrCat["E_scat_onsphere_cart_N",Sprintf["%g.pos",pe]]],
                Name StrCat["E_scat_onsphere_cart_N",Sprintf["%g",pe]]];
        EndIf
        For po In {1:p_max}
          Print[feN~{pe}~{po}[SurfInt] , OnGlobal, Format Table, File StrCat[myDir,StrCat[StrCat["feN_",Sprintf["pe%g",pe]],Sprintf["po%g.dat",po]]]];
          Print[fhN~{pe}~{po}[SurfInt] , OnGlobal, Format Table, File StrCat[myDir,StrCat[StrCat["fhN_",Sprintf["pe%g",pe]],Sprintf["po%g.dat",po]]]];
        EndFor      
      }
    }
  EndFor
}

DefineConstant[
  R_ = {"res_VPWall_helmholtz_vector", Name "GetDP/1ResolutionChoices", Visible 1},
  C_ = {"-solve -pos -petsc_prealloc 200 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps", Name "GetDP/9ComputeCommand", Visible 1},
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}];

