import h5py
import sys
import os
import numpy as np
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import utilities.tmatrix_tools as tools
from config.config import Config
import meep as mp
import treams
from t_computation.scat import get_start_geometry
def h5save (c:Config, t, file=None,  name = "t_data", description = ""):
    speed_of_light = 299792458.
    basis = treams.SphericalWaveBasis.default(c.l_max)
    pols = np.array(tools.translate_pols(basis.pol, "parity"))
    if not mp.am_really_master():
        return
    keywords = c.keywords
    geometry_resolution = np.array([*c.start_geometry.shape]) / c.object_size
    n_pad = ((c.no_pml_size - c.object_size) * geometry_resolution / 2.0).astype(
        int
    )
    geo_dat = get_start_geometry(c)
    eps_grid = np.asarray(geo_dat.eps_grid, dtype=float)    
    no_pad_mask = np.asarray(geo_dat.no_pad_mask, dtype=float)
    with h5py.File(c.path_output + f"{name}.tmat.h5", "w") as f:
        tools.base_data(
            f,
            tmatrices=t,
            name=name,
            description=description,
            keywords=keywords,
            freqs=c.frequencies*speed_of_light,
            ftype="frequency",
            funit="MHz",
            modes=(basis.l,basis.m, pols),
            format_version="v1"
            )     
        f.create_group("embedding")
        for index, (name, description, keywords, epsilon, mu) in enumerate(zip(
            c.material_names,
            c.material_descriptions,
            c.material_keywords,
            np.append(c.eps_embedding, c.material),
            c.mu
        )):
            if index == 0:
                tools.isotropic_material(
                f[f"embedding"],
                epsilon = c.eps_embedding,
                name=c.material_names[0],
                description=c.material_descriptions[0],
                keywords=c.material_keywords[0],
                #epsilon=c.start_geometry,
                mu=c.mu[0]
            )                
            else:
                if  isinstance(c.material, list):
                    scatname =  f"scatterer_{index}"
                    f.create_group(f"{scatname}/geometry/")
                    tools.geometry_shape(
                    f[f"{scatname}/geometry"],
                    shape = c.shape[index-1],
                    params = c.params[index-1],
                    lunit="um",
                    )
                    f[f"{scatname}/geometry/"].attrs["unit"] = "um"
                    f[f"{scatname}/geometry/position"] = c.positions[index-1]
                    f[f"{scatname}/geometry/position"].attrs["unit"] = "um"
                    
                else:
                    scatname = f"scatterer"
                    # f.create_group(f"{scatname}/geometry/")
                    f[f"{scatname}/geometry/mesh"] = h5py.SoftLink("/computation/mesh")
                    f[f"{scatname}/geometry/"].attrs["unit"] = "um"
                    tools.geometry_shape(
                    f[f"{scatname}/geometry"],
                    shape = c.shape,
                    params = c.params,
                    lunit="um",
                    )
                f.create_group(f"{scatname}/material")
                tools.isotropic_material(
                    f[f"{scatname}/material"],
                    name=c.material_names[index],
                    description=c.material_descriptions[index],
                    keywords=c.material_keywords[index],
                    epsilon=c.material,
                    mu=c.mu[index]
                )
        f.create_group("computation/files")
        f.create_group("computation/mesh")
        f["computation/mesh"].create_dataset('eps_grid', data=eps_grid, compression='gzip')
        f["computation/mesh"].attrs["description"] = "eps_grid is the voxel grid of permittivity distribution inside a box of sidelength defined by the parameter domain in the method_parameters. This does not include the PML."

        tools.computation_data(
            f["computation"],
            name="MEEP",
            description="Using plane wave illuminations",
            keywords="",
            keys= { "resolution": c.resolution,
                    "domain": c.object_size,
                    "sim_amount_mult": c.sim_amount_mult, 
                    "dpml": c.dpml,
                    "f_src": c.f_src, 
                    "cpu_cores_per_simulation": c.cpu_cores_per_simulation
                },
            method="Finite Difference Time Domain, FDTD", 
            software=f"MEEP={mp.__version__}, python={sys.version_info[0]}.{sys.version_info[1]}.{sys.version_info[2]}, numpy={np.__version__}, h5py={h5py.__version__}",            
        )
        with open(__file__, "r") as scriptfile:
            f[f"computation/files/{os.path.basename(__file__)}"] = scriptfile.read()
        f[f"computation/files/source_script.py"] = file
        f["mesh"] = h5py.SoftLink("/computation/mesh")
