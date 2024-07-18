import numpy as np
import treams.io
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import utilities.tmatrix_tools as tools
import h5py
import onelab


###Parameters
################################################################################
lambdas = np.linspace(750, 1250, 101)  # steps of 5 nm
shape = "cylinder"
radius = 250.0
height = 300.0
params = {"radius": radius, "height": height}
eps_cyl = 2.5**2
eps_emb = 1.0
material_names = ["Air, Vacuum", "TiO2, Titanium Dioxide"]
material_descriptions = ["", ""]
material_keywords = ["nondispersive", "nondispersive, dielectric"]
material_values = [eps_emb, eps_cyl]
mu = [1.0, 1.0]
lmax = 4
ms = 12.0
pmax = lmax * lmax + 2 * lmax
name = f"Tmatrix_cylinder_TiO2_lm_{lmax}_ms_{int(ms)}"

###Functions onelab > treams
################################################################################


def load_getdp_integral(filename):
    return np.loadtxt(filename)[1] + 1j * np.loadtxt(filename)[2]


def _get_postprocess(wl):
    myDir = "resT_wl_" + str(wl) + "/"
    Tmatrix = np.zeros((2 * pmax, 2 * pmax), dtype=complex)
    for ne in range(1, lmax + 1):
        for no in range(1, lmax + 1):
            for me in range(-ne, ne + 1):
                for mo in range(-no, no + 1):
                    pe = int(ne * (ne + 1) - me)
                    po = int(no * (no + 1) - mo)
                    ke = int(ne * (ne + 1) + me)
                    ko = int(no * (no + 1) + mo)
                    Tmatrix[2 * ke - 2, 2 * ko - 2] = load_getdp_integral(
                        myDir + "feN_pe%gpo%g.dat" % (pe, po)
                    )
                    Tmatrix[2 * ke - 1, 2 * ko - 1] = load_getdp_integral(
                        myDir + "fhM_pe%gpo%g.dat" % (pe, po)
                    )
                    Tmatrix[2 * ke - 2, 2 * ko - 1] = load_getdp_integral(
                        myDir + "feM_pe%gpo%g.dat" % (pe, po)
                    )
                    Tmatrix[2 * ke - 1, 2 * ko - 2] = load_getdp_integral(
                        myDir + "fhN_pe%gpo%g.dat" % (pe, po)
                    )
    return Tmatrix


###Create and initialize onelab client
################################################################################

# create a new onelab client
c = onelab.client(__file__)

# get Gmsh and GetDP locations from Gmsh options

mygmsh = c.getString("General.ExecutableFileName")
mygetdp = ""

for s in range(9):
    n = c.getString("Solver.Name" + str(s))
    if n == "GetDP":
        mygetdp = c.getString("Solver.Executable" + str(s))
        break

if not len(mygetdp):
    c.sendError("This appears to be the first time you are trying to run GetDP")
    c.sendError(
        "Please run a GetDP model interactively once with Gmsh to "
        + "initialize the solver location"
    )
    exit(0)

c.sendInfo("Will use gmsh={0} and getdp={1}".format(mygmsh, mygetdp))

# create a onelab variable for the model name
mymodel = c.defineString("Model name", value="scattererTmatrix")

# we're done if we don't do the actual calculation
if c.action == "check":
    exit(0)

# get model file names with correct path
mymodel_geo = c.getPath(mymodel + ".geo")
mymodel_msh = c.getPath(mymodel + ".msh")
mymodel_pro = c.getPath(mymodel + ".pro")

c.sendInfo("The path to access my model is {0}".format(mymodel_geo))


###Perform the job
################################################################################
c.setNumber("1Geometry/00Scatterer shape", value=2)
c.setNumber("1Geometry/01cylinder X-radius [nm]", value=radius)
c.setNumber("1Geometry/02cylinder Y-radius [nm]", value=radius)
c.setNumber("1Geometry/03cylinder height [nm]", value=height)
c.setNumber(
    "3Electromagnetic parameters/00scatterer permittivity (real) []",
    value=np.real(eps_cyl),
)
c.setNumber(
    "3Electromagnetic parameters/01scatterer permittivity (imag) []",
    value=np.imag(eps_cyl),
)
c.setNumber("3Electromagnetic parameters/08n_max integer", value=lmax)
c.setNumber("4Mesh size and PMLs parameters/02mesh size", value=ms)

tmats = []
for lambda0 in lambdas:
    c.setNumber("3Electromagnetic parameters/04wavelength [nm]", value=int(lambda0))

    # run gmsh as a subclient
    c.runSubClient("myGmsh", mygmsh + " " + mymodel_geo + " -3 -o " + mymodel_msh)

    # run getdp as a subclient
    c.runSubClient(
        "myGetDP",
        mygetdp
        + " "
        + mymodel_pro
        + " -pre res_VPWall_helmholtz_vector -msh "
        + mymodel_msh
        + " -cal -petsc_prealloc 200",
    )

    # build T-matrix
    single = treams.PhysicsArray(
        _get_postprocess(int(lambda0)),
        basis=treams.SphericalWaveBasis.default(lmax),
        k0=2 * np.pi / lambda0,
        material=1,
        modetype=("singular", "regular"),
        poltype="parity",
    )
    print(single)
    tmats.append(single)


# Store data in the standard format .tmat.h5
with h5py.File("cylinder.tmat.h5", "w") as fobj:
    pol = tmats[0].basis.pol
    pol = np.where(pol == 1, "electric", "magnetic")
    tools.base_data(
        fobj,
        tmatrices=tmats,
        name="TiO2 cylinder in air",
        description=(
            f"Reference TiO2 cylinder from the tmatrix_data_format example with {radius} nm radius and {height} nm height."
        ),
        keywords="czinfinity,mirrorxyz,passive, reciprocal, lossless",
        freqs=lambdas,
        ftype="vacuum_wavelength",
        funit=tools.LENGTHS[1e-9],
        modes=(tmats[0].basis.l, tmats[0].basis.m, pol),
        format_version="v1",
    )
    for index, (name, description, keywords, epsilon, mu) in enumerate(
        zip(
            material_names,
            material_descriptions,
            material_keywords,
            material_values,
            mu,
        )
    ):
        if index == 0:
            fobj.create_group(f"embedding")
            tools.isotropic_material(
                fobj[f"embedding"],
                name=name,
                description=description,
                keywords=keywords,
                epsilon=epsilon,
                mu=mu,
            )
        else:
            if len(material_values) == 2:
                scatname = "scatterer"
            else:
                scatname = f"scatterer_{index}"
            fobj.create_group(f"{scatname}")
            fobj.create_group(f"{scatname}/material")
            tools.isotropic_material(
                fobj[f"{scatname}/material"],
                name=name,
                description=description,
                keywords=keywords,
                epsilon=epsilon,
                mu=mu,
            )
    fobj.create_group("computation/files")
    tools.computation_data(
        fobj["computation"],
        name="ONELAB",
        description="This computation used adapted version of https://gitlab.onelab.info/doc/models/-/tree/master/ElectromagneticScattering, which was detailed in: Demésy, Guillaume, Jean-Claude Auger, and Brian Stout. 'Scattering matrix of arbitrarily shaped objects: combining finite elements and vector partial waves.' JOSA A 35.8 (2018): 1401-1409.",
        method="FEM, Finite Element Method",
        keys={"mesh_size": ms},
        software=f"onelab=1.0, python={sys.version.split()[0]}, numpy={np.__version__}",
        # meshfile=os.path.join(working_dir, "grid.jcm"),
        lunit=tools.LENGTHS[1e-9],
    )
    with open(__file__, "r") as scriptfile:
        fobj[f"computation/files/{os.path.basename(__file__)}"] = scriptfile.read()

    with open("scattererTmatrix.msh", "r") as mf:
        mesh = mf.read()

    fobj["computation/mesh.msh"] = mesh
    fobj.create_group(f"{scatname}/geometry")
    fobj[f"{scatname}/geometry/mesh.msh"] = h5py.SoftLink("/computation/mesh.msh")

    fobj[f"{scatname}/geometry"].attrs["unit"] = tools.LENGTHS[1e-9]
    tools.geometry_shape(
        fobj[f"{scatname}/geometry"],
        shape=shape,
        params=params,
        lunit=tools.LENGTHS[1e-9],
    )
    fobj["mesh"] = h5py.SoftLink("/computation/mesh.msh")
