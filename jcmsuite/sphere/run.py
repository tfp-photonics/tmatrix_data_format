"""Template for a T-matrix calculation with JCMsuite at the example of a sphere."""
# -> set your JCMROOT installation directory
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import h5py
import jcmwave
import meshgen as meshparser
import numpy as np
import tmatrix_tools


def run():
    # Setup
    jcmt_pattern = "sphere"
    working_dir = "tmp"

    keys = {
        "uol": 1e-9,  # unit of length (1 = m, 1e-9 = nm)
        "epsilon": [1, 4],  # First entry embedding
        "mu": [1, 1],
        "radius": 100,
        "degree_max": 5,
    }
    keys_method = {
        "domain_radius": 125,
        "domain_z": 250,
        "maximum_sidelength_domain": 10,
        "maximum_sidelength_object": 5,
        "precision": 1e-7,
        "max_refinements": 3,
        "fem_degree": 2,
    }
    keys.update(keys_method)
    uof = 1e12  # unit of frequency (1 = Hz, 1e12 = THz)
    freqs = np.linspace(400, 750, 10)
    c0 = 299792458 / (uof * keys["uol"])  # speed of light
    print("freqs", freqs)
    material_names = ["Vacuum", "Custom"]
    material_descriptions = ["", ""]
    material_keywords = ["non-dispersive", "non-dispersive"]

    jcmwave.geo(".", keys, jcmt_pattern=jcmt_pattern, working_dir=working_dir)
    jcmwave.daemon.shutdown()  # Make sure no daemon is running
    jcmwave.daemon.add_workstation(Multiplicity=8)

    jobids = []
    for i, freq in enumerate(freqs):
        keys["lambda0"] = c0 * keys["uol"] / freqs[i]
        jobid = jcmwave.solve(
            "project.jcmp",
            keys=keys,
            working_dir=f"{working_dir}/job{i}",
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
    with h5py.File(jcmt_pattern + ".h5", "w") as fobj:
        tmatrix_tools.base_data(
            fobj,
            tmatrices=tmats,
            name="Sphere",
            description="A simple example of a single sphere at 10 frequencies.",
            keywords="czinfinity,mirrorxyz,passive, reciprocal",
            freqs=freqs,
            ftype="frequency",
            funit=tmatrix_tools.FREQUENCIES[uof],
            modes=(ls, ms, pols),
            modes_inc=modes_inc,
            format_version="draft",
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
            description="Using the built-in post-process 'MultipoleExpansion'",
            method="FEM, Finite Element Method",
            program_version="jcmsuite=" + jcmwave.__private.version,
            keys=keys_method,
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
            shape="sphere",
            params=keys,
            meshfile=mesh,
            lunit=tmatrix_tools.LENGTHS[keys["uol"]],
            working_dir=working_dir,
        )
        fobj["mesh.msh"] = h5py.SoftLink("/computation/mesh.msh")


if __name__ == "__main__":
    run()
