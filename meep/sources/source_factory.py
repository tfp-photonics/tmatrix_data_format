import meep as mp
import numpy as np

from config.config import Config
from sources.source_data import SourceData
from utility.vector_geometry import index_to_meep_e_field_component

def get_plane_wave_source(c: Config, k_src_dir = mp.Vector3(0.0,0.0,1.0), pol_src_dir = mp.Vector3(1.0,0.0,0.0)) -> SourceData:
    if(np.abs(k_src_dir.dot(pol_src_dir)) > 0.0001 ):
        raise ValueError("Wave vector and polarization vector of a plane wave must be orthogonal.")
    k_src_unit_vec = k_src_dir/k_src_dir.norm()
    pol_src_unit_vec = pol_src_dir/pol_src_dir.norm()
    k_src_abs = 2 * np.pi * c.f_src  # normal units and c=1 so k=2*pi*f
    k_src = k_src_abs*k_src_unit_vec

    def pw_amp(k, x0):
        def _pw_amp(x):
            return np.exp(1j * k.dot(x + x0))

        return _pw_amp

    source_plane_normal_axis = np.argmax(np.abs(k_src))
    source_plane_normal_direction = -np.sign(k_src[source_plane_normal_axis])
    source_size = np.copy(c.with_pml_size)
    source_size[source_plane_normal_axis] = 0
    source_center = np.zeros(3)
    source_center[source_plane_normal_axis] = (
        source_plane_normal_direction * c.no_pml_size[source_plane_normal_axis] / 2.0
    )
    sources = []
    for pol_component_index, pol_component in enumerate(pol_src_unit_vec):
        source = mp.Source(
            mp.GaussianSource(
                frequency=c.f_src,
                fwidth=c.df_src,
                is_integrated=True,
            ),
            center=mp.Vector3(*source_center),
            size=mp.Vector3(*source_size),
            component=index_to_meep_e_field_component(pol_component_index),
            amplitude=pol_component,
            amp_func=pw_amp(k_src, mp.Vector3()),
        )
        sources.append(source)

    src_data = SourceData(
        meep_sources=sources, k_src=np.array([*k_src]), pol_src=np.array([*pol_src_unit_vec])
    )
    return src_data
