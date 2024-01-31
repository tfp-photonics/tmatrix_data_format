import h5py
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(os.path.dirname(current_dir))

sys.path.append(parent_dir)
print(sys.path)
import tmatrix_tools
# import tmat_data_format.tmatrix_data_format_repo.tmatrix_tools as tmatrix_tools
from config.config import Config
import meep as mp
import treams

def persist_T(c:Config, t, name = "t_data", description = ""):
    basis = treams.SphericalWaveBasis.default(c.l_max)
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
            freqs=c.frequencies,
            ftype="frequency",
            funit="MHz",
            modes=(basis.l,basis.m,basis.pol)
            )
        f.create_group("materials/voxel_grid")
        tmatrix_tools.isotropic_material(
                f["materials/voxel_grid"],
                name="Relative permittivity distribution",
                description=f"Voxel grid of relative permittivity values of the scatterer.",
                keywords=keywords,
                epsilon=c.start_geometry,
                mu=1,
                embedding=c.eps_embedding,
        )
        f.create_group("geometry")
        tmatrix_tools.geometry_shape(
            f["geometry"],
            shape="rectangular_cuboid",
            params={
                "lengthx":c.object_size[0],
                "lengthy":c.object_size[1],
                "lengthz":c.object_size[2],
                "object_size": c.object_size,
                "dpml": c.dpml,
                "resolution": c.resolution,
                "simulation_amount_mult": c.sim_amount_mult,
                "f_src": c.f_src,
                "df_src": c.df_src
            },
            lunit="micrometers",
        )