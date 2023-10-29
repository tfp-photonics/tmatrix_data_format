/*
 * sphere_simple.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

public class sphere_simple {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.param().set("l_max", "5", "Maximum degree");
    model.param().set("l_in", "1", "Degree incident wave");
    model.param().set("m_in", "0", "Order incident wave");
    model.param().set("p_in", "1", "Polarization incident wave");
    model.param().set("l_out", "1", "Degree scattered wave");
    model.param().group().create("par2");
    model.param("par2").set("r_decomp", "150[nm]", "Radius of the decomposition");
    model.param("par2").set("r_domain", "200[nm]", "Radius of the domain");
    model.param("par2").set("d_pml", "1000[nm]", "Thickness PML");
    model.param("par2").set("r_obj", "100[nm]", "Radius of the object");
    model.param("default").label("Parameters Decomposition");
    model.param("par2").label("Parameters Geometry");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 2);

    model.result().table().create("tbl1", "Table");

    model.func().create("an1", "Analytic");
    model.func().create("an2", "Analytic");
    model.func().create("an3", "Analytic");
    model.func().create("an4", "Analytic");
    model.func().create("an5", "Analytic");
    model.func().create("an6", "Analytic");
    model.func().create("an7", "Analytic");
    model.func().create("an8", "Analytic");
    model.func().create("an9", "Analytic");
    model.func().create("an10", "Analytic");
    model.func().create("an11", "Analytic");
    model.func().create("an12", "Analytic");
    model.func().create("an13", "Analytic");
    model.func().create("an14", "Analytic");
    model.func().create("an15", "Analytic");
    model.func().create("an16", "Analytic");
    model.func().create("an17", "Analytic");
    model.func().create("an18", "Analytic");
    model.func().create("an19", "Analytic");
    model.func().create("an20", "Analytic");
    model.func().create("an21", "Analytic");
    // model.func().create("an22", "Analytic");
    // model.func().create("an23", "Analytic");
    // model.func().create("an24", "Analytic");
    // model.func().create("an25", "Analytic");
    // model.func().create("an26", "Analytic");
    // model.func().create("an27", "Analytic");
    // model.func().create("an28", "Analytic");
    model.func("an1").label("Spherical Bessel First Kind");
    model.func("an1").set("funcname", "sphbesselj");
    model.func("an1").set("expr", "sqrt(pi / (2 * x)) * besselj(l + 0.5, x)");
    model.func("an1").set("args", new String[] {"l", "x"});
    model.func("an2").label("Spherical Hankel Second Kind");
    model.func("an2").set("funcname", "sphhankelh2");
    model.func("an2").set("expr", "sqrt(pi / (2 * x)) * (besselj(l + 0.5, x) - i * bessely(l + 0.5, x))");
    model.func("an2").set("args", new String[] {"l", "x"});
    model.func("an2").set("complex", true);
    model.func("an3").label("Psi function");
    model.func("an3").set("funcname", "psi2");
    model.func("an3").set("expr", "(l + 1) * sphhankelh2(l, x) / x - sphhankelh2(l + 1, x)");
    model.func("an3").set("args", new String[] {"l", "x"});
    model.func("an3").set("complex", true);
    model.func("an4").label("Pi function");
    model.func("an4").set("funcname", "pilm");
    model.func("an4").set("expr", "-0.5 * (legendre(l + 1, m + 1, cos(theta)) + (l - m + 1) * (l - m + 2) * legendre(l + 1, m - 1, cos(theta)))");
    model.func("an4").set("args", new String[] {"l", "m", "theta"});
    model.func("an5").label("Tau function");
    model.func("an5").set("funcname", "taulm");
    model.func("an5").set("expr", "(legendre(l, m + 1, cos(theta)) - (l + m) * (l - m + 1) * legendre(l, m - 1, cos(theta))) * 0.5");
    model.func("an5").set("args", new String[] {"l", "m", "theta"});
    model.func("an6").label("Prefactor VSW");
    model.func("an6").set("funcname", "gammalm");
    model.func("an6").set("expr", "i * sqrt((2 * l + 1) * factorial(l - m) / (4 * pi * l * (l + 1) * factorial(l + m)))");
    // model.func("an6").set("expr", "i * sqrt((2 * l + 1) * gamma(l - m + 1) / (4 * pi * l * (l + 1) * gamma(l + m + 1)))");
    model.func("an6").set("args", new String[] {"l", "m"});
    model.func("an6").set("complex", true);
    model.func("an7").label("Vector Spherical Xrho");
    model.func("an7").set("funcname", "Xlmrho");
    model.func("an7").set("expr", "gammalm(l, m) * exp(i * m * phi) * i * cos(theta) * pilm(l, m, theta)");
    model.func("an7").set("args", new String[] {"l", "m", "theta", "phi"});
    model.func("an7").set("complex", true);
    model.func("an8").label("Vector Spherical Xphi");
    model.func("an8").set("funcname", "Xlmphi");
    model.func("an8").set("expr", "-gammalm(l, m) * exp(i * m * phi) * taulm(l, m, theta)");
    model.func("an8").set("args", new String[] {"l", "m", "theta", "phi"});
    model.func("an8").set("complex", true);
    model.func("an9").label("Vector Spherical Xz");
    model.func("an9").set("funcname", "Xlmz");
    model.func("an9").set("expr", "-gammalm(l, m) * exp(i * m * phi) * i * m * legendre(l, m, cos(theta))");
    model.func("an9").set("args", new String[] {"l", "m", "theta", "phi"});
    model.func("an9").set("complex", true);
    model.func("an10").label("Vector Spherical r cross Xrho");
    model.func("an10").set("funcname", "rXlmrho");
    model.func("an10").set("expr", "gammalm(l, m) * exp(i * m * phi) * cos(theta) * taulm(l, m, theta)");
    model.func("an10").set("args", new String[] {"l", "m", "theta", "phi"});
    model.func("an10").set("complex", true);
    model.func("an11").label("Vector Spherical r cross Xphi");
    model.func("an11").set("funcname", "rXlmphi");
    model.func("an11").set("expr", "gammalm(l, m) * exp(i * m * phi) * i * pilm(l, m, theta)");
    model.func("an11").set("args", new String[] {"l", "m", "theta", "phi"});
    model.func("an11").set("complex", true);
    model.func("an12").label("Vector Spherical r cross Xz");
    model.func("an12").set("funcname", "rXlmz");
    model.func("an12").set("expr", "-gammalm(l, m) * exp(i * m * phi) * sin(theta) * taulm(l, m, theta)");
    model.func("an12").set("args", new String[] {"l", "m", "theta", "phi"});
    model.func("an12").set("complex", true);
    model.func("an13").label("Vector Spherical RgMrho");
    model.func("an13").set("funcname", "RgMrho");
    model.func("an13").set("expr", "sphbesselj(l, x) * Xlmrho(l, m, theta, phi)");
    model.func("an13").set("args", new String[] {"l", "m", "x", "theta", "phi"});
    model.func("an13").set("complex", true);
    model.func("an14").label("Vector Spherical RgMphi");
    model.func("an14").set("funcname", "RgMphi");
    model.func("an14").set("expr", "sphbesselj(l, x) * Xlmphi(l, m, theta, phi)");
    model.func("an14").set("args", new String[] {"l", "m", "x", "theta", "phi"});
    model.func("an14").set("complex", true);
    model.func("an15").label("Vector Spherical RgMz");
    model.func("an15").set("funcname", "RgMz");
    model.func("an15").set("expr", "sphbesselj(l, x) * Xlmz(l, m, theta, phi)");
    model.func("an15").set("args", new String[] {"l", "m", "x", "theta", "phi"});
    model.func("an15").set("complex", true);
    model.func("an16").label("Vector Spherical RgNrho");
    model.func("an16").set("funcname", "RgNrho");
    model.func("an16").set("expr", "rXlmrho(l, m, theta, phi) * ((l + 1) * sphbesselj(l, x) / x - sphbesselj(l + 1, x)) + gammalm(l, m) * exp(i * m * phi) * sin(theta) * l * (l + 1) * legendre(l, m, cos(theta)) * sphbesselj(l, x) / x");
    model.func("an16").set("args", new String[] {"l", "m", "x", "theta", "phi"});
    model.func("an16").set("complex", true);
    model.func("an17").label("Vector Spherical RgNphi");
    model.func("an17").set("funcname", "RgNphi");
    model.func("an17").set("expr", "rXlmphi(l, m, theta, phi) * ((l + 1) * sphbesselj(l, x) / x - sphbesselj(l + 1, x))");
    model.func("an17").set("args", new String[] {"l", "m", "x", "theta", "phi"});
    model.func("an17").set("complex", true);
    model.func("an18").label("Vector Spherical RgNz");
    model.func("an18").set("funcname", "RgNz");
    model.func("an18").set("expr", "rXlmz(l, m, theta, phi) * ((l + 1) * sphbesselj(l, x) / x- sphbesselj(l + 1, x)) + gammalm(l, m) * exp(i * m * phi) * cos(theta) * l * (l + 1) * legendre(l, m, cos(theta)) * sphbesselj(l, x) / x");
    model.func("an18").set("args", new String[] {"l", "m", "x", "theta", "phi"});
    model.func("an18").set("complex", true);
    model.func("an19").label("Vector Spherical RgArho");
    model.func("an19").set("funcname", "RgArho");
    model.func("an19").set("expr", "sqrt(0.5) * (RgNrho(l, m, x, theta, phi) + p * RgMrho(l, m, x, theta, phi))");
    model.func("an19").set("args", new String[] {"l", "m", "p", "x", "theta", "phi"});
    model.func("an19").set("complex", true);
    model.func("an20").label("Vector Spherical RgAphi");
    model.func("an20").set("funcname", "RgAphi");
    model.func("an20").set("expr", "sqrt(0.5) * (RgNphi(l, m, x, theta, phi) + p * RgMphi(l, m, x, theta, phi))");
    model.func("an20").set("args", new String[] {"l", "m", "p", "x", "theta", "phi"});
    model.func("an20").set("complex", true);
    model.func("an21").label("Vector Spherical RgAz");
    model.func("an21").set("funcname", "RgAz");
    model.func("an21").set("expr", "sqrt(0.5) * (RgNz(l, m, x, theta, phi) + p * RgMz(l, m, x, theta, phi))");
    model.func("an21").set("args", new String[] {"l", "m", "p", "x", "theta", "phi"});
    model.func("an21").set("complex", true);
    // model.func("an22").label("Legendre 1");
    // model.func("an22").set("funcname", "legendre1");
    // model.func("an22").set("expr", "if(m == 0, x, -sqrt(1 - x^2))");
    // model.func("an22").set("args", new String[] {"m", "x"});
    // model.func("an23").label("Legendre 2");
    // model.func("an23").set("funcname", "legendre2");
    // model.func("an23").set("expr", "if(m == 0, (1.5 * x^2 - 0.5), if(m == 1, -3 * x * sqrt(1 - x^2), 3 * (1 - x^2)))");
    // model.func("an23").set("args", new String[] {"m", "x"});
    // model.func("an24").label("Legendre 3");
    // model.func("an24").set("funcname", "legendre3");
    // model.func("an24").set("expr", "if(m == 0, 2.5 * x^3 - 1.5 * x, if(m == 1, -(7.5 * x^2 - 1.5) * sqrt(1 - x^2), if(m == 2, 15 * x * (1 - x^2), -15 * sqrt(1 - x^2)^3)))");
    // model.func("an24").set("args", new String[] {"m", "x"});
    // model.func("an25").label("Legendre 4");
    // model.func("an25").set("funcname", "legendre4");
    // model.func("an25").set("expr", "if(m == 0, 4.375 * x^4 - 3.75 * x^2 + 0.375, if(m == 1, -(17.5 * x^3 - 7.5 * x) * sqrt(1 - x^2), if(m == 2, (52.5 * x^2 - 7.5) * (1 - x^2), if(m == 3, -105 * x * sqrt(1 - x^2)^3, 105 * (1 - x^2)^2))))");
    // model.func("an25").set("args", new String[] {"m", "x"});
    // model.func("an26").label("Legendre 5");
    // model.func("an26").set("funcname", "legendre5");
    // model.func("an26").set("expr", "if(m == 0, 7.875 * x^5 - 8.75 * x^3 + 1.875 * x, if(m == 1, -(39.375 * x^4 - 26.25 * x^2 + 1.875) * sqrt(1 - x^2), if(m == 2, (157.5 * x^3 - 52.5 * x) * (1 - x^2), if(m == 3, -(472.5 * x^2 - 52.5) * sqrt(1 - x^2)^3, if(m == 4, 945 * x * (1 - x^2)^2, -945 * sqrt(1 - x^2)^5)))))");
    // model.func("an26").set("args", new String[] {"m", "x"});
    // model.func("an27").label("Legendre positive");
    // model.func("an27").set("funcname", "legendrepos");
    // model.func("an27").set("expr", "if(m > l, 0, if(l == 0, 1, if(l == 1, legendre1(m, x), if(l == 2, legendre2(m, x), if(l == 3, legendre3(m, x), if(l == 4, legendre4(m, x), if(l == 5, legendre5(m, x), -99999999)))))))");
    // model.func("an27").set("args", new String[] {"l", "m", "x"});
    // model.func("an28").label("Associated Legendre");
    // model.func("an28").set("funcname", "legendre");
    // model.func("an28").set("expr", "if(m < 0, (-1)^m * gamma(l + m + 1) * legendrepos(l, abs(m), x) / gamma(l - m + 1), legendrepos(l, m, x))");
    // model.func("an28").set("args", new String[] {"l", "m", "x"});

    model.component("comp1").geom("geom1").axisymmetric(true);

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").create("c1", "Circle");
    model.component("comp1").geom("geom1").feature("c1").label("Domain");
    model.component("comp1").geom("geom1").feature("c1").set("rot", -90);
    model.component("comp1").geom("geom1").feature("c1").set("layername", new String[] {"Layer 1", "Layer 2"});
    model.component("comp1").geom("geom1").feature("c1").set("layer", new String[] {"d_pml", "r_domain - r_decomp"});
    model.component("comp1").geom("geom1").feature("c1").set("r", "r_domain + d_pml");
    model.component("comp1").geom("geom1").feature("c1").set("angle", 180);
    model.component("comp1").geom("geom1").create("c2", "Circle");
    model.component("comp1").geom("geom1").feature("c2").label("Object");
    model.component("comp1").geom("geom1").feature("c2").set("rot", -90);
    model.component("comp1").geom("geom1").feature("c2").set("r", "r_obj");
    model.component("comp1").geom("geom1").feature("c2").set("angle", 180);
    model.component("comp1").geom("geom1").run();
    model.component("comp1").geom("geom1").run("fin");

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("rs", "sqrt(r * r + z * z)", "Radius");
    model.component("comp1").variable("var1").set("theta", "atan2(r, z)", "Polar angle");
    model.component("comp1").variable("var1").set("n_domain", "sqrt(mat1.def.epsilonr_iso * mat1.def.mur_iso)", "Refractive index domain");
    model.component("comp1").variable("var1").set("kappa_domain", "mat1.pg1.kappa", "Kappa domain");
    model.component("comp1").variable().create("var2");
    model.component("comp1").variable("var2").set("xp", "emw.k0 * r_decomp * (n_domain + kappa_domain)");
    model.component("comp1").variable("var2").set("xm", "emw.k0 * r_decomp * (n_domain - kappa_domain)");
    model.component("comp1").variable("var2").set("intX", "2 * pi * intop1(sin(theta) * (Xlmrho(l_out, m_in, theta, 0) * withsol('sol2', emw.relEr, setval(l_in, l_in), setval(m_in, m_in), setval(p_in, p_in), setval(freq, freq)) + Xlmphi(l_out, m_in, theta, 0) * withsol('sol2', emw.relEphi, setval(l_in, l_in), setval(m_in, m_in), setval(p_in, p_in), setval(freq, freq)) + Xlmz(l_out, m_in, theta, 0) * withsol('sol2', emw.relEz, setval(l_in, l_in), setval(m_in, m_in), setval(p_in, p_in), setval(freq, freq)))) / r_decomp");
    model.component("comp1").variable("var2").set("intrX", "2 * pi * intop1(sin(theta) * (rXlmrho(l_out, m_in, theta, 0) * withsol('sol2', emw.relEr, setval(l_in, l_in), setval(m_in, m_in), setval(p_in, p_in), setval(freq, freq)) + rXlmphi(l_out, m_in, theta, 0) * withsol('sol2', emw.relEphi, setval(l_in, l_in), setval(m_in, m_in), setval(p_in, p_in), setval(freq, freq)) + rXlmz(l_out, m_in, theta, 0) * withsol('sol2', emw.relEz, setval(l_in, l_in), setval(m_in, m_in), setval(p_in, p_in), setval(freq, freq)))) / r_decomp");
    model.component("comp1").variable("var2").set("shh2p", "sphhankelh2(l_out, xp)");
    model.component("comp1").variable("var2").set("shh2m", "sphhankelh2(l_out, xm)");
    model.component("comp1").variable("var2").set("psi2p", "psi2(l_out, xp)");
    model.component("comp1").variable("var2").set("psi2m", "psi2(l_out, xm)");
    model.component("comp1").variable("var2").set("determ", "shh2p * psi2m + shh2m * psi2p");
    model.component("comp1").variable("var2").set("ap", "sqrt(2) * conj((psi2m * intX + shh2m * intrX) / determ) * 1[m/V]");
    model.component("comp1").variable("var2").set("am", "sqrt(2) * conj((-psi2p * intX + shh2p * intrX) / determ) * 1[m/V]");

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat1").propertyGroup().create("pg1", "Biisotropy");
    model.component("comp1").material("mat2").selection().set(4);
    model.component("comp1").material("mat2").propertyGroup().create("pg1", "Biisotropy");

    model.component("comp1").cpl().create("intop1", "Integration");
    model.component("comp1").cpl("intop1").selection().geom("geom1", 1);
    model.component("comp1").cpl("intop1").selection().set(13, 16);

    model.component("comp1").coordSystem().create("pml1", "PML");
    model.component("comp1").coordSystem("pml1").selection().set(1, 6);

    model.component("comp1").physics().create("emw", "ElectromagneticWaves", "geom1");

    model.component("comp1").variable("var1").label("Variables Simulation");
    model.component("comp1").variable("var2").label("Variables Decomposition");

    model.component("comp1").material("mat1").label("Material Domain");
    model.component("comp1").material("mat1").propertyGroup("def")
        .set("relpermittivity", new String[] {"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat1").propertyGroup("def").set("relpermittivity_symmetry", "0");
    model.component("comp1").material("mat1").propertyGroup("def")
        .set("relpermeability", new String[] {"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat1").propertyGroup("def").set("relpermeability_symmetry", "0");
    model.component("comp1").material("mat1").propertyGroup("def")
        .set("electricconductivity", new String[] {"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat1").propertyGroup("def").set("electricconductivity_symmetry", "0");
    model.component("comp1").material("mat1").propertyGroup("pg1").set("kappa", "0");
    model.component("comp1").material("mat1").propertyGroup("pg1").descr("kappa", "Chirality parameter");
    model.component("comp1").material("mat1").propertyGroup("pg1").set("chi", "0");
    model.component("comp1").material("mat1").propertyGroup("pg1").descr("chi", "Non-reciprocity parameter");
    model.component("comp1").material("mat2").label("Material Object");
    model.component("comp1").material("mat2").propertyGroup("def")
        .set("relpermittivity", new String[] {"4", "0", "0", "0", "4", "0", "0", "0", "4"});
    model.component("comp1").material("mat2").propertyGroup("def").set("relpermittivity_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("def")
        .set("relpermeability", new String[] {"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat2").propertyGroup("def").set("relpermeability_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("def")
        .set("electricconductivity", new String[] {"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    model.component("comp1").material("mat2").propertyGroup("def").set("electricconductivity_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("pg1").set("kappa", "0");
    model.component("comp1").material("mat2").propertyGroup("pg1").descr("kappa", "Chirality parameter");
    model.component("comp1").material("mat2").propertyGroup("pg1").set("chi", "0");
    model.component("comp1").material("mat2").propertyGroup("pg1").descr("chi", "Non-reciprocity parameter");

    model.component("comp1").physics("emw").prop("BackgroundField").set("SolveFor", "scatteredField");
    model.component("comp1").physics("emw").prop("BackgroundField")
        .set("Eb", new String[][] {
            {"conj(RgArho(l_in, m_in, p_in, emw.k0 * conj(n_domain + p_in * kappa_domain) * rs, theta, 0)) * 1[V/m]"},
            {"conj(RgAphi(l_in, m_in, p_in, emw.k0 * conj(n_domain + p_in * kappa_domain) * rs, theta, 0)) * 1[V/m]"},
            {"conj(RgAz(l_in, m_in, p_in, emw.k0 * conj(n_domain + p_in * kappa_domain) * rs, theta, 0)) * 1[V/m]"}});
    model.component("comp1").physics("emw").prop("outofplanewavenumber").set("mFloquet", "m_in");
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.dHdtr", new String[] {"(emw.murinvrr * (emw.dBdtr - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Er / c_const) + emw.murinvrphi * (emw.dBdtphi - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Ephi / c_const) + emw.murinvrz * (emw.dBdtz - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Ez / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.dHdtphi", new String[] {"(emw.murinvphir * (emw.dBdtr - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Er / c_const) + emw.murinvphiphi * (emw.dBdtphi - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Ephi / c_const) + emw.murinvphiz * (emw.dBdtz - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Ez / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.dHdtz", new String[] {"(emw.murinvzr * (emw.dBdtr - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Er / c_const) + emw.murinvzphi * (emw.dBdtphi - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Ephi / c_const) + emw.murinvzz * (emw.dBdtz - (material.pg1.chi + i * material.pg1.kappa) * emw.iomega * emw.Ez / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.Dr", new String[] {"epsilon0_const * emw.Er + emw.Pr + (material.pg1.chi - i * material.pg1.kappa) * emw.Hr / c_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.Dphi", new String[] {"epsilon0_const * emw.Ephi + emw.Pphi + (material.pg1.chi - i * material.pg1.kappa) * emw.Hphi / c_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.Dz", new String[] {"epsilon0_const * emw.Ez + emw.Pz + (material.pg1.chi - i * material.pg1.kappa) * emw.Hz / c_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.Hr", new String[] {"(emw.murinvrr * (emw.Br - (material.pg1.chi + i * material.pg1.kappa) * emw.Er / c_const) + emw.murinvrphi * (emw.Bphi - (material.pg1.chi + i * material.pg1.kappa) * emw.Ephi / c_const) + emw.murinvrz * (emw.Bz - (material.pg1.chi + i * material.pg1.kappa) * emw.Ez / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.Hphi", new String[] {"(emw.murinvphir * (emw.Br - (material.pg1.chi + i * material.pg1.kappa) * emw.Er / c_const) + emw.murinvphiphi * (emw.Bphi - (material.pg1.chi + i * material.pg1.kappa) * emw.Ephi / c_const) + emw.murinvphiz * (emw.Bz - (material.pg1.chi + i * material.pg1.kappa) * emw.Ez / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").feature("wee1").featureInfo("info")
        .set("emw.Hz", new String[] {"(emw.murinvzr * (emw.Br - (material.pg1.chi + i * material.pg1.kappa) * emw.Er / c_const) + emw.murinvzphi * (emw.Bphi - (material.pg1.chi + i * material.pg1.kappa) * emw.Ephi / c_const) + emw.murinvzz * (emw.Bz - (material.pg1.chi + i * material.pg1.kappa) * emw.Ez / c_const)) / mu0_const"});

    model.component("comp1").physics("emw").featureInfo("info")
        .set("emw.Dbr", new String[] {"epsilon0_const * emw.Ebr + emw.Pbr + (material.pg1.chi - i * material.pg1.kappa) * emw.Hbr / c_const"});
    model.component("comp1").physics("emw").featureInfo("info")
        .set("emw.Dbphi", new String[] {"epsilon0_const * emw.Ebphi + emw.Pbphi + (material.pg1.chi - i * material.pg1.kappa) * emw.Hbphi / c_const"});
    model.component("comp1").physics("emw").featureInfo("info")
        .set("emw.Dbz", new String[] {"epsilon0_const * emw.Ebz + emw.Pbz + (material.pg1.chi - i * material.pg1.kappa) * emw.Hbz / c_const"});
    model.component("comp1").physics("emw").featureInfo("info")
        .set("emw.Hbr", new String[] {"(emw.murinvrr * (emw.Bbr - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebr / c_const) + emw.murinvrphi * (emw.Bbphi - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebphi / c_const) + emw.murinvrz * (emw.Bbz - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebz / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").featureInfo("info")
        .set("emw.Hbphi", new String[] {"(emw.murinvphir * (emw.Bbr - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebr / c_const) + emw.murinvphiphi * (emw.Bbphi - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebphi / c_const) + emw.murinvphiz * (emw.Bbz - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebz / c_const)) / mu0_const"});
    model.component("comp1").physics("emw").featureInfo("info")
        .set("emw.Hbz", new String[] {"(emw.murinvzr * (emw.Bbr - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebr / c_const) + emw.murinvzphi * (emw.Bbphi - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebphi / c_const) + emw.murinvzz * (emw.Bbz - (material.pg1.chi + i * material.pg1.kappa) * emw.Ebz / c_const)) / mu0_const"});

    model.study().create("std1");
    model.study("std1").create("param3", "Parametric");
    model.study("std1").create("param2", "Parametric");
    model.study("std1").create("param", "Parametric");
    model.study("std1").create("freq", "Frequency");
    model.study("std1").setGenPlots(false);
    model.study().create("std2");
    model.study("std2").create("param4", "Parametric");
    model.study("std2").create("param3", "Parametric");
    model.study("std2").create("param2", "Parametric");
    model.study("std2").create("param", "Parametric");
    model.study("std2").create("freq", "Frequency");
    model.study("std2").feature("freq").set("activate", new String[] {"emw", "off"});

    model.study("std1").feature("param").label("Parametric Sweep l_in");
    model.study("std1").feature("param").set("pname", new String[] {"l_in"});
    model.study("std1").feature("param").set("plistarr", new String[] {"range(1, 1, l_max)"});
    model.study("std1").feature("param").set("punit", new String[] {""});
    model.study("std1").feature("param2").label("Parametric Sweep m_in");
    model.study("std1").feature("param2").set("pname", new String[] {"m_in"});
    model.study("std1").feature("param2").set("plistarr", new String[] {"range(-l_in, 1, l_in)"});
    model.study("std1").feature("param2").set("punit", new String[] {""});
    model.study("std1").feature("param3").label("Parametric Sweep p_in");
    model.study("std1").feature("param3").set("pname", new String[] {"p_in"});
    model.study("std1").feature("param3").set("plistarr", new String[] {"1, -1"});
    model.study("std1").feature("param3").set("punit", new String[] {""});
    model.study("std1").feature("freq").set("punit", "THz");
    model.study("std1").feature("freq").set("plist", "range(400, 1, 750)");
    model.study("std1").feature("freq").set("plot", false);
    model.study("std1").feature("freq").set("ngen", 5);
    model.study("std2").feature("param").label("Parametric Sweep l_out");
    model.study("std2").feature("param").set("pname", new String[] {"l_out"});
    model.study("std2").feature("param").set("plistarr", new String[] {"range(1, 1, l_max)"});
    model.study("std2").feature("param").set("punit", new String[] {""});
    model.study("std2").feature("param2").label("Parametric Sweep l_in");
    model.study("std2").feature("param2").set("pname", new String[] {"l_in"});
    model.study("std2").feature("param2").set("plistarr", new String[] {"range(1, 1, l_max)"});
    model.study("std2").feature("param2").set("punit", new String[] {""});
    model.study("std2").feature("param3").label("Parametric Sweep m_in");
    model.study("std2").feature("param3").set("pname", new String[] {"m_in"});
    model.study("std2").feature("param3").set("plistarr", new String[] {"range(-min(l_in, l_out), 1, min(l_in, l_out))"});
    model.study("std2").feature("param3").set("punit", new String[] {""});
    model.study("std2").feature("param4").label("Parametric Sweep p_in");
    model.study("std2").feature("param4").set("pname", new String[] {"p_in"});
    model.study("std2").feature("param4").set("plistarr", new String[] {"1, -1"});
    model.study("std2").feature("param4").set("punit", new String[] {""});
    model.study("std2").feature("freq").set("punit", "THz");
    model.study("std2").feature("freq").set("plist", "range(400, 1, 750)");
    model.study("std2").feature("freq").set("preusesol", "no");

    model.result().table("tbl1").set("storetable", "onfile");
    model.result().table("tbl1").set("filename", "tmatrix_coeffs.txt");
    model.result().numerical().create("gev1", "EvalGlobal");
    model.result().numerical("gev1").set("expr", new String[] {"ap", "am"});
    model.result().numerical("gev1").set("table", "tbl1");


    // Comment this to avoid exporting the mesh
    model.component("comp1").mesh("mesh1").run();
    model.mesh("mesh1").export("mesh1.mphtxt");

    // Comment this to not start the run immediately
    model.study("std1").run();
    model.study("std2").run();
    model.result().numerical("gev1").set("data", "dset4");
    model.result().numerical("gev1").setResult();

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
