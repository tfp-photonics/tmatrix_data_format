"""Template for a T-matrix calculation with JCMsuite at the example of a cylinder."""
# -> set your JCMROOT installation directory
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
jcm_root = '/scratch/local/nasadova/JCMsuite' # -> set your JCMROOT installation directory
sys.path.append(os.path.join(jcm_root, 'ThirdPartySupport', 'Python'))
import h5py
import jcmwave
import numpy as np
import meshgen as meshparser
import tmatrix_tools

def run():
    # Setup
    jcmt_pattern = "cylinder"
    working_dir = "tmp"

    keys = {
        "uol": 1e-9,  # unit of length (1 = m, 1e-9 = nm)
        "epsilon": [1, 6.25],
        "mu": [1, 1],
        "radius": 250,
        "height": 300,
        "degree_max": 4
        }
    
    keys_method = {     
        "domain_radius": 275,
        "domain_z": 350,           
        "maximum_sidelength_domain": 15,
        "maximum_sidelength_object": 6,
        "precision": 1e-7,
        "max_refinements": 3,
        "fem_degree": 2,
    }
    keys.update(keys_method)
    uof = 1e12  # unit of frequency (1 = Hz, 1e12 = THz)
    freqs = np.linspace(240, 400, 101)
    c0 = 299792458 / (uof * keys["uol"])  # speed of light
    material_names = ["Vacuum, Air", "TiO2, Titanium dioxide"]
    material_descriptions = ["", "A constant refractive index of 2.5 is used."]
    material_keywords = ["non-dispersive", "non-dispersive"]

    jcmwave.geo(".", keys, jcmt_pattern=jcmt_pattern, working_dir=working_dir)
    jcmwave.daemon.shutdown()  # Make sure no daemon is running
    jcmwave.daemon.add_workstation(Multiplicity=20) #How many instances to run in parallel (each uses one license)

    jobids = []
    for i, freq in enumerate(freqs):
        keys["lambda0"] = c0 * keys["uol"] / freq
        jobid = jcmwave.solve(
            "project.jcmp",
            keys=keys,
            working_dir=os.path.join(working_dir, f"job{i}"),
            jcmt_pattern=jcmt_pattern,
        )
        jobids.append(jobid)
    results, logs = jcmwave.daemon.wait(jobids)
    jcmwave.daemon.shutdown()

    tmats, _, ls, ms, pols = tmatrix_tools.extract_tmatrix_jcm(results)
    
    ls_inc, ms_inc, pols_inc = tmatrix_tools.jcm_defaultmodes(keys["degree_max"])
    tmats, pols, pols_inc = tmatrix_tools.jcm_pol_reorder(tmats, pols, pols_inc, keys["degree_max"])
    if np.all(ls == ls_inc) and np.all(ms == ms_inc) and np.all(pols == pols_inc):
        modes_inc = None
    else:
        modes_inc = ls, ms, pols
    with h5py.File(jcmt_pattern + "_tio2.tmat.h5", "w") as fobj:
        tmatrix_tools.base_data(
            fobj,
            tmatrices=tmats,
            name="Cylinder",
            description=(
                f"A simple cylinder with {keys['radius']} nm radius"
                f" and {keys['height']} nm height."
            ),
            keywords="czinfinity,mirrorxyz,passive, reciprocal",
            freqs=freqs,
            ftype="frequency",
            funit=tmatrix_tools.FREQUENCIES[uof],
            modes=(ls, ms, pols),
            modes_inc=modes_inc,
            format_version="v1",
        )
        for index, (name, description, keywords, epsilon, mu) in enumerate(zip(
            material_names,
            material_descriptions,
            material_keywords,
            keys["epsilon"],
            keys["mu"]
        )):
            if index == 0:
                fobj.create_group(f"embedding")
                tmatrix_tools.isotropic_material(
                    fobj[f"embedding"],
                    name=name,
                    description=description,
                    keywords=keywords,
                    epsilon=epsilon,
                    mu=mu
                )
            else:
                if len(keys["epsilon"]) == 2:
                    scatname = "scatterer"
                else:
                    scatname = f"scatterer_{index}"
                fobj.create_group(f"{scatname}")
                fobj.create_group(f"{scatname}/material")
                tmatrix_tools.isotropic_material(
                    fobj[f"{scatname}/material"],
                    name=name,
                    description=description,
                    keywords=keywords,
                    epsilon=epsilon,
                    mu=mu
                )
        fobj.create_group("computation/files")
        tmatrix_tools.computation_data(
            fobj["computation"],
            name="JCMsuite",
            description="Using the built-in post-process 'MultipoleExpansion'. Rotationally symmetric 3D problem is computed in 2D",
            method="FEM, Finite Element Method",
            keys=keys_method,
            program_version=f"jcmsuite={jcmwave.__private.version}, python={sys.version.split()[0]}, numpy={np.__version__}, h5py={h5py.__version__}",
            meshfile=os.path.join(working_dir, "grid.jcm"),
            lunit=tmatrix_tools.LENGTHS[keys["uol"]],
        )
        tmatrix_tools.jcm_files(
            fobj["computation"],
            jcmt_pattern,
            ".",
        )
        with open(__file__, "r") as scriptfile:
            fobj[f"computation/files/{os.path.basename(__file__)}"] = scriptfile.read()
        with open(meshparser.__file__, "r") as scriptfile:
            fobj[f"computation/files/{os.path.basename(meshparser.__file__)}"] = scriptfile.read()
        fobj.create_group(f"{scatname}/geometry")
        mesh = meshparser.Mesh(
            {
                domain: domain.tag
                for domain in meshparser.read(
                    os.path.join(working_dir, "grid.jcm")
                ).domains
                if domain.tag > 1 and domain.dim == 2 # change to 3 for 3D
            }
        )
        fobj[f"{scatname}/geometry"].attrs["unit"] = tmatrix_tools.LENGTHS[keys["uol"]]
        tmatrix_tools.geometry_shape(
            fobj[f"{scatname}/geometry"],
            shape="cylinder",
            params=keys,
            meshfile=mesh,
            lunit=tmatrix_tools.LENGTHS[keys["uol"]],
            working_dir=working_dir,
        )
        fobj["mesh"] = h5py.SoftLink("/computation/mesh.msh")



if __name__ == "__main__":
    run()
