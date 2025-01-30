import io
import os
import re
import sys

import h5py
import numpy as np

try:
    import meshgen as meshparser
except ImportError as _err:
    meshparser = _err


LENGTHS = {
    1e-24: "ym",
    1e-21: "zm",
    1e-18: "am",
    1e-15: "fm",
    1e-12: "pm",
    1e-9: "nm",
    1e-6: "um",
    1e-3: "mm",
    1e-2: "cm",
    1e-1: "dm",
    1: "m",
    1e1: "dam",
    1e2: "hm",
    1e3: "km",
    1e6: "Mm",
    1e9: "Gm",
    1e12: "Tm",
    1e15: "Pm",
    1e18: "Em",
    1e21: "Zm",
    1e24: "Ym",
}

INVLENGTHS = {
    1e24: r"ym^{-1}",
    1e21: r"zm^{-1}",
    1e18: r"am^{-1}",
    1e15: r"fm^{-1}",
    1e12: r"pm^{-1}",
    1e9: r"nm^{-1}",
    1e6: r"um^{-1}",
    1e3: r"mm^{-1}",
    1e2: r"cm^{-1}",
    1e1: r"dm^{-1}",
    1: r"m^{-1}",
    1e-1: r"dam^{-1}",
    1e-2: r"hm^{-1}",
    1e-3: r"km^{-1}",
    1e-6: r"Mm^{-1}",
    1e-9: r"Gm^{-1}",
    1e-12: r"Tm^{-1}",
    1e-15: r"Pm^{-1}",
    1e-18: r"Em^{-1}",
    1e-21: r"Zm^{-1}",
    1e-24: r"Ym^{-1}",
}

FREQUENCIES = {
    1e-24: "yHz",
    1e-21: "zHz",
    1e-18: "aHz",
    1e-15: "fHz",
    1e-12: "pHz",
    1e-9: "nHz",
    1e-6: "uHz",
    1e-3: "mHz",
    1e-2: "cHz",
    1e-1: "dHz",
    1: "Hz",
    1e1: "daHz",
    1e2: "hHz",
    1e3: "kHz",
    1e6: "MHz",
    1e9: "GHz",
    1e12: "THz",
    1e15: "PHz",
    1e18: "EHz",
    1e21: "ZHz",
    1e24: "YHz",
}


def translate_pols(pols, poltype):
    """
    This assumes particular ordering of initial ones and zeros,
    does not ensure correct final order.
    """
    if poltype == "parity":
        dct = {
            "magnetic": 0,
            "te": 0,
            "M": 0,
            "electric": 1,
            "tm": 1,
            "N": 1,
            0: "magnetic",
            1: "electric",
        }
    elif poltype == "helicity":
        dct = {
            "negative": 0,
            "minus": 0,
            "positive": 1,
            "plus": 1,
            -1: "negative",
            0: "negative",
            1: "positive",
        }
    return [dct[pol] for pol in pols]


def _name_descr_kw(fobj, name, description="", keywords=""):
    for key, val in [
        ("name", name),
        ("description", description),
        ("keywords", keywords),
    ]:
        val = str(val)
        if val != "":
            fobj.attrs[key] = val


def base_data(
    fobj,
    *,
    tmatrices,
    name,
    description="",
    keywords="",
    freqs,
    ftype="frequency",
    funit="Hz",
    modes,
    modes_inc=None,
    format_version="v0.0.1",
):
    fobj["tmatrix"] = np.asarray(tmatrices)

    _name_descr_kw(fobj, name, description, keywords)

    if ftype not in (
        "frequency",
        "angular_frequency",
        "vacuum_wavelength",
        "vacuum_wavenumber",
        "angular_vacuum_wavenumber",
    ):
        raise ValueError(f"invalid frequency/wavenumber/wavelength type {ftype}")
    fobj[ftype] = np.asarray(freqs)
    fobj[ftype].attrs["unit"] = funit

    if modes_inc is None:
        fobj["modes/l"] = np.asarray(modes[0])
        fobj["modes/m"] = np.asarray(modes[1])
        fobj.create_dataset(
            "modes/polarization", data=[*map(str, modes[2])], dtype=h5py.string_dtype()
        )
    else:
        fobj["modes/l_scattered"] = np.asarray(modes[0])
        fobj["modes/m_scattered"] = np.asarray(modes[1])
        fobj.create_dataset(
            "modes/polarization_scattered",
            data=[*map(str, modes[2])],
            dtype=h5py.string_dtype(),
        )
        fobj["modes/l_incident"] = np.asarray(modes_inc[0])
        fobj["modes/m_incident"] = np.asarray(modes_inc[1])
        fobj.create_dataset(
            "modes/polarization_incident",
            data=[*map(str, modes_inc[2])],
            dtype=h5py.string_dtype(),
        )

    fobj.attrs["storage_format_version"] = format_version


def _check_name(fobj, name, exact=True):
    if (exact and fobj.name != name) or not fobj.name.startswith(name):
        raise ValueError(f"require name '{name}'")


def isotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    reference="",
    interpolation="",
    index=None,
    impedance=None,
    epsilon=None,
    mu=None,
):
    # _check_name(fobj.parent, "/scatterer/material")
    if reference != "":
        fobj.attrs["reference"] = reference
    if interpolation != "":
        fobj.attrs["interpolation"] = interpolation
    _name_descr_kw(fobj, name, description, keywords)
    if index is impedance is epsilon is mu is None:
        TypeError("missing at least one of: 'index', 'impedance', 'epsilon', 'mu'")
    if not (index is impedance is None) and not (epsilon is mu is None):
        TypeError("'epsilon' and 'mu' exclude the use of 'index' and 'impedance'")
    for param, val in [
        ("relative_permittivity", epsilon),
        ("relative_permeability", mu),
        ("refractive_index", index),
        ("relative_impedance", impedance),
    ]:
        if val is not None:
            fobj[param] = np.asarray(val)


def anisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    index=None,
    epsilon=None,
    mu=None,
    reference="",
    interpolation="",
    index_inner_dims=2,
    epsilon_inner_dims=2,
    mu_inner_dims=2,
    coordinates="Cartesian",
):
    isotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        reference=reference,
        interpolation=interpolation,
        index=index,
        epsilon=epsilon,
        mu=mu,
    )
    for param, val, inner_dims in [
        ("relative_permittivity", epsilon, epsilon_inner_dims),
        ("relative_permeability", mu, mu_inner_dims),
        ("refractive_index", index, index_inner_dims),
    ]:
        if val is not None:
            fobj[param].attrs["inner_dims"] = int(inner_dims)
            fobj[param].attrs["coordinate_system"] = str(coordinates)


def biisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    reference="",
    interpolation="",
    epsilon=None,
    mu=None,
    kappa=None,
    chi=None
):
    isotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        reference=reference,
        interpolation=interpolation,
        epsilon=epsilon,
        mu=mu,
    )
    for param, val in [
        ("chirality", kappa),
        ("nonreciprocity", chi),
    ]:
        if val is not None:
            fobj[param] = np.asarray(val)


def bianisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    reference="",
    interpolation="",
    epsilon=None,
    mu=None,
    kappa=None,
    chi=None,
    epsilon_inner_dims=2,
    mu_inner_dims=2,
    kappa_inner_dims=2,
    chi_inner_dims=2,
    coordinates="Cartesian",
):
    biisotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        reference=reference,
        interpolation=interpolation,
        epsilon=epsilon,
        mu=mu,
        kappa=kappa,
        chi=chi,
    )
    for param, val, inner_dims in [
        ("relative_permittivity", epsilon, epsilon_inner_dims),
        ("relative_permeability", mu, mu_inner_dims),
        ("chirality", kappa, kappa_inner_dims),
        ("nonreciprocity", chi, mu_inner_dims),
    ]:
        if val is not None:
            fobj[param].attrs["inner_dims"] = int(inner_dims)
            fobj[param].attrs["coordinate_system"] = str(coordinates)


def bianisotropic_material_from_tensor(
    fobj,
    *,
    name,
    description="",
    keywords="",
    reference="",
    interpolation="",
    bianisotropy,
    inner_dims=2,
    coordinates="Cartesian",
):
    _check_name(fobj.parent, "/materials")
    if reference != "":
        fobj.attrs["reference"] = reference
    if interpolation != "":
        fobj.attrs["interpolation"] = interpolation
    _name_descr_kw(fobj, name, description, keywords)
    fobj["bianisotropy"] = np.asarray(bianisotropy)
    fobj["bianisotropy"].attrs["inner_dims"] = int(inner_dims)
    fobj["bianisotropy"].attrs["coordinate_system"] = str(coordinates)


def mesh_data(fobj, meshfile, lunit="m"):
    if isinstance(meshparser, ImportError):
        return
    if isinstance(meshfile, str):
        mesh = meshparser.Mesh.read(meshfile)
    elif isinstance(meshfile, meshparser.Mesh):
        mesh = meshfile
    else:
        raise ValueError(f"invalid type of meshfile: {type(meshfile)}")
    buffer = io.StringIO()
    mesh.write_gmsh(buffer)
    buffer.seek(0)
    fobj["mesh.msh"] = buffer.read()
    fobj["mesh.msh"].attrs["unit"] = lunit


def computation_data(
    fobj,
    *,
    name="",
    description="",
    keywords="",
    method="",
    software="",
    keys="",
    meshfile=None,
    lunit="m",
    working_dir=".",
):
    _check_name(fobj, "/computation")
    _name_descr_kw(fobj, name, description, keywords)
    fobj.attrs["software"] = software
    if "method" != "":
        fobj.attrs["method"] = method
    fobj.create_group("method_parameters")
    for key in keys:
        fobj["method_parameters"][f"{key}"] = np.asarray(keys[key])
    if meshfile is not None:
        mesh_data(fobj, meshfile, lunit)
        with open(os.path.join(working_dir, meshfile), "rb") as fobj_mesh:
            fobj[f"files/{os.path.split(meshfile)[1]}"] = np.void(fobj_mesh.read())
    with open(__file__, "r") as scriptfile:
        fobj[f"files/{os.path.basename(__file__)}"] = scriptfile.read()


GEOMETRY_PARAMS = {
    "sphere": ("radius",),
    "ellipsoid": ("radiusx", "radiusy", "radiusz"),
    "superellipsoid": ("radiusx", "radiusy", "radiusz", "e_parm", "n_parm"),
    "spheroid": ("radiusxy", "radiusz"),
    "cylinder": ("radius", "height"),
    "cone": ("radius_top", "radius_bottom", "height"),
    "torus": ("radius_major", "radius_minor"),
    "cube": ("length",),
    "rectangular_cuboid": ("lengthx", "lengthy", "lengthz"),
    "helix": (
        "radius_helix",
        "radius_wire",
        "number_turns",
        "pitch",
        "handedness",
        "termination",
    ),
    "ring": (
        "radius_major", 
        "radius_minor", 
        "height"
    ),
    "convex_polyhedron": ("points",), 
    "pyramid": (
        "n_edges", 
        "radius", 
        "height", 
        "angle", 
        "apex_shift"
        ),
    "regular_prism": (
        "n_edges", 
        "radius", 
        "height", 
        "shift"
    ),
    "wedge": (
        "lengthx", 
        "lengthy", 
        "lengthz", 
        "deltax",
        "deltay"
    )
}


def geometry_shape(
    fobj,
    *,
    shape,
    params,
    name="",
    description="",
    keywords="",
    meshfile=None,
    lunit="m",
    working_dir=".",
):

    _name_descr_kw(fobj, name, description, keywords)
    for param in GEOMETRY_PARAMS[shape]:
        fobj[param] = params[param]
    fobj.attrs["shape"] = str(shape)
    if meshfile is not None:
        mesh_data(fobj, meshfile, lunit)
        if isinstance(meshfile, str):
            with open(os.path.join(working_dir, meshfile), "rb") as fobj_mesh:
                fobj[os.path.split(meshfile)[1]] = np.void(fobj_mesh.read())
