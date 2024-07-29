"""Template for a T-matrix calculation with COMSOL at the example of a sphere."""

import os
import re
import h5py
import numpy as np
import sys
sys.path.append('../..')
sys.path.append('../../..')
import tmatrix_tools as tools
def run():
    working_dir = "."
    comsolname = "sphere_simple"
    material_names = ["Embedding", "Sphere"]
    material_descriptions = ["", ""]
    material_keywords = ["", ""]
    epsilon = [1, 4]
    mu = [1, 1]
    kappa = [0., 0.]
    (
        tmats,
        _,
        ls,
        ms,
        pols,
        ls_inc,
        ms_inc,
        pols_inc,
        params,
    ) = tools.extract_tmatrix_comsol(os.path.join(working_dir, "tmatrix_coeffs.txt"))

    if np.all(ls == ls_inc) and np.all(ms == ms_inc) and np.all(pols == pols_inc):
        modes_inc = None
    else:
        modes_inc = ls_inc, ms_inc, pols_inc
    for key, val in params.items():
        freq = re.match(r"^freq \(([a-zA-Z]*Hz)\)$", key)
        if freq is not None:
            funit = freq[1]
            freqs = np.asarray(val, float)
            break
    else:
        raise ValueError("frequency definition not found")
    with h5py.File(comsolname + ".h5", "w") as fobj:
        tools.base_data(
            fobj,
            tmatrices=tmats,
            name="Sphere",
            description=f"A simple example of a single sphere at {len(freqs)} frequencies.",
            keywords="czinfinity,mirrorxyz, passive, reciprocal",
            freqs=freqs,
            ftype="frequency",
            funit=funit,
            modes=(ls, ms, pols),
            modes_inc=modes_inc,
            format_version="draft",
        )

        for index, (name, description, keywords, epsilon, mu, kappa) in enumerate(zip(
            material_names,
            material_descriptions,
            material_keywords,
            epsilon,
            mu, 
            kappa
        )):
            if index == 0:
                fobj.create_group(f"embedding")
                tools.isotropic_material(
                    fobj[f"embedding"],
                    name=name,
                    description=description,
                    keywords=keywords,
                    epsilon=epsilon,
                    mu=mu, 
                    kappa=kappa
                )
            else:
                if len(epsilon) == 2:
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
                    mu=mu
                )
        fobj.create_group("computation/files")
        fobj.create_group(f"{scatname}/geometry")
        fobj[f"{scatname}/geometry"].attrs["unit"] = "nm"
        tools.geometry_shape(
            fobj[f"{scatname}/geometry"],
            shape="sphere",
            params="radius",
            meshfile=None,
            lunit="nm",
            working_dir=working_dir,
        )
        fobj.create_group("computation/files")
        tools.computation_data(
            fobj["computation"],
            name="COMSOL",
            description="Using a custom multipole decomposition.",
            keywords="FEM, Finite Element Method",
            program_version="comsol=6.1.0.522",
            meshfile=None,
            # meshfile="mesh1.mphtxt",
            
        )
        with open(os.path.join(working_dir, f"{comsolname}.java"), "r") as javafile:
            fobj[f"computation/files/{comsolname}.java"] = javafile.read()
        with open(__file__, "r") as scriptfile:
            fobj[f"computation/files/{os.path.basename(__file__)}"] = scriptfile.read()


if __name__ == "__main__":
    run()
