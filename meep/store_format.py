import h5py
import sys
import os
import numpy as np
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.append(parent_dir)
import tmatrix_tools
# import tmat_data_format.tmatrix_data_format_repo.tmatrix_tools as tmatrix_tools
from config.config import Config
import meep as mp
import treams

def persist_T(c:Config, t, name = "t_data", description = ""):
    speed_of_light = 299792458.
    basis = treams.SphericalWaveBasis.default(c.l_max)
    pols = np.array(tmatrix_tools.translate_pols(basis.pol, "parity"))
    if not mp.am_really_master():
        return
    keywords = "MEEP, FDTD, Voxel"
    with h5py.File(c.path + f"{name}.h5", "w") as f:
        tmatrix_tools.base_data(
            f,
            tmatrices=t,
            name=name,
            description=description,
            keywords=keywords,
            freqs=c.frequencies*speed_of_light,
            ftype="frequency",
            funit="MHz",
            modes=(basis.l,basis.m, pols)
            )

        f.create_group("geometry")
        if isinstance(c.shape, list):
            for i in range(len(c.shape)):
                print('i', i)
                f.create_group(f"geometry/object{i}")
                tmatrix_tools.geometry_shape(
                f["geometry"][f"object{i}"],
                shape = c.shape[i],
                params = c.params[i],
                lunit="um",
                )
                f.create_group(f"materials/object{i}")  
                tmatrix_tools.isotropic_material(
                f['materials'][f'object{i}'],
                epsilon = c.material[i],
                name = "",
                description = "",
                #name="Relative permittivity distribution",
                #description=f"Voxel grid of relative permittivity values of the scatterer.",
                keywords=keywords,
                #epsilon=c.start_geometry,
                mu=1,
                #embedding=c.eps_embedding,
            
            )
            f["geometry"]["positions"] = c.positions
            f['materials'].parent.parent["embedding"] = h5py.SoftLink(f['materials'].name)

        else:
            f.create_group(f"geometry/object")
            tmatrix_tools.geometry_shape(
            f["geometry"],
            shape = c.shape,
            params = c.params,

            lunit="um",
        )
            f.create_group("materials/object")
            tmatrix_tools.isotropic_material(
                #f["materials/voxel_grid"],
                f['materials/object'],
                epsilon = c.material,
                name = "",
                description = "",
                #name="Relative permittivity distribution",
                #description=f"Voxel grid of relative permittivity values of the scatterer.",
                keywords=keywords,
                #epsilon=c.start_geometry,
                mu=1,
                embedding=c.eps_embedding,
            
            )
        f.create_group("computation")
        tmatrix_tools.computation_data(
            f["computation"],
            name="MEEP",
            description="Using plane wave illumination",
            keywords="FDTD",
                                                                                          )
