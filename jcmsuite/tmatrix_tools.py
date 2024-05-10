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

def jcm_defaultmodes(lmax):
    res = np.array(
        [
            [l, m, pol]
            for l in range(1, lmax + 1)  # noqa: E741
            for m in range(-l, l + 1)
            for pol in range(2)
        ]
    ).T
    return res[0], res[1], translate_pols(res[2], "parity")

def jcm_pol_reorder(tmats, pols, pols_inc, lmax):
    """
    Default jcm order is magnetic/electric in parity basis 
    Modify to reverse order

    Args:
        lmax (int): maximal degree
        tmats (array): T-matrices
        pols (array): polarization array of scattered modes
        pols_inc (array): polarization array of incident modes


    Returns:
        Tmatrices and polarizations of incident and scattered modes        
    """
    if (pols[0] == "magnetic"):
        reorder  = [
                    2 * (l * l + l - 1 + m) + 1 - p
                    for l in range(1, lmax + 1)
                    for m in range(-l, l + 1)
                    for p in range(2)
                ]
        if tmats.ndim == 2:
            tmats = np.expand_dims(tmats, axis=0)
        for i, t in enumerate(tmats):
            ti = t[:, reorder]
            tmats[i][:, :] = ti[reorder, :]
    pols_new = np.zeros_like(pols)
    pols_inc_new = np.zeros_like(pols_inc)
    pols_new[::2] = pols[1::2]
    pols_new[1::2] = pols[::2]
    pols_inc_new[1::2] = pols_inc[::2]
    pols_inc_new[::2] = pols_inc[1::2]
    return tmats, pols_new, pols_inc_new

def translate_pols(pols, poltype):
    if poltype == "parity":
        dct = {
            "magnetic": 0,
            "te": 0,
            "M": 0,
            "electric": 1,
            "tm": 1,
            "N": 0,
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


def _get_postprocess(result, name):
    processes = [process for process in result if process.get("title", "") == name]
    if len(processes) != 1:
        raise ValueError(f"can't find post-process '{name}'")
    return processes[0]


def extract_tmatrix_jcm(results):
    """Turn result list of JCM to list of T-matrices.

    Args:
        results (list): JCM result.
        lmax (int): Maximal degree.
        npostprocess (int, optional): Index of the MultipoleExpansion postprocess.

    Returns:
        A list of the T-matrices
    """
    nresults = len(results)
    if nresults == 0:
        return np.zeros((0, 0, 0), complex), np.zeros(0)
    angular_frequencies = np.empty(nresults)
    angular_frequencies[:] = np.nan
    expansion = _get_postprocess(results[0], "ElectricFieldStrength_MultipoleExpansion")
    ls = np.asarray(expansion["n"], int)
    ms = np.asarray(expansion["m"], int)
    pols = np.asarray(expansion["Type"], int)
    keys_inc = expansion["ElectricFieldStrength"].keys()
    tmatrices = np.empty((nresults, len(ls), len(keys_inc)), complex)
    for i, result in enumerate(results):
        expansion = _get_postprocess(result, "ElectricFieldStrength_MultipoleExpansion")
        for j, col in expansion["ElectricFieldStrength"].items():
            tmatrices[i, :, j] = col
        angular_frequencies[i] = np.real(expansion["header"]["FDOmega"])
        if (
            np.any(ls != expansion["n"])
            or np.any(ms != expansion["m"])
            or np.any(pols != expansion["Type"])
            or keys_inc != expansion["ElectricFieldStrength"].keys()
        ):
            raise ValueError("mixed post-process types for 'MultipoleExpansion'")
    return (
        tmatrices,
        angular_frequencies,
        ls,
        ms,
        np.array(translate_pols(pols, "parity")),
    )


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
    fobj.attrs["created_with"] = "python=" + sys.version.split()[0]


def _check_name(fobj, name, exact=True):
    if (exact and fobj.name != name) or not fobj.name.startswith(name):
        raise ValueError(f"require name '{name}'")


def isotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    index=None,
    impedance=None,
    epsilon=None,
    mu=None,
    embedding=False,
):
    # _check_name(fobj.parent, "/scatterer/material")

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
    if embedding:
        fobj.parent.parent["embedding"] = h5py.SoftLink(fobj.name)


def anisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    index=None,
    epsilon=None,
    mu=None,
    index_inner_dims=2,
    epsilon_inner_dims=2,
    mu_inner_dims=2,
    coordinates="Cartesian",
    embedding=False,
):
    isotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        index=index,
        epsilon=epsilon,
        mu=mu,
        embedding=embedding,
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
    epsilon=None,
    mu=None,
    kappa=None,
    chi=None,
    embedding=False,
):
    isotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        epsilon=epsilon,
        mu=mu,
        embedding=embedding,
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
    epsilon=None,
    mu=None,
    kappa=None,
    chi=None,
    epsilon_inner_dims=2,
    mu_inner_dims=2,
    kappa_inner_dims=2,
    chi_inner_dims=2,
    coordinates="Cartesian",
    embedding=False,
):
    biisotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        epsilon=epsilon,
        mu=mu,
        kappa=kappa,
        chi=chi,
        embedding=embedding,
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
    bianisotropy,
    inner_dims=2,
    coordinates="Cartesian",
    embedding=False,
):
    _check_name(fobj.parent, "/materials")

    _name_descr_kw(fobj, name, description, keywords)
    fobj["bianisotropy"] = np.asarray(bianisotropy)
    fobj["bianisotropy"].attrs["inner_dims"] = int(inner_dims)
    fobj["bianisotropy"].attrs["coordinate_system"] = str(coordinates)
    if embedding:
        fobj.parent.parent["embedding"] = h5py.SoftLink(fobj.name)


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
    program_version="",
    keys="",
    meshfile=None,
    lunit="m",
    working_dir=".",
):
    _check_name(fobj, "/computation")
    _name_descr_kw(fobj, name, description, keywords)
    fobj.attrs["software"] = program_version
    if "method"!="":
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

def jcm_files(fobj, jcmt_pattern=None, working_dir="."):
    _check_name(fobj, "/computation")
    jcm_file(fobj["files"], "project", jcmt_pattern, "jcmp", working_dir)
    for filename in ("boundary_condition", "layout", "materials", "sources"):
        jcm_file(fobj["files"], filename, jcmt_pattern, "jcm", working_dir)


def jcm_file(fobj, name, pattern=None, ending="jcm", working_dir="."):
    _check_name(fobj, "/computation/files")
    for parts in [[name, pattern, ending], [name, ending]]:
        fullname = ".".join([i for i in parts if i is not None])
        if os.path.isfile(os.path.join(working_dir, fullname + "t")):
            with open(os.path.join(working_dir, fullname + "t"), "r") as jcmfile:
                fobj[fullname + "t"] = jcmfile.read()
            break
        elif os.path.isfile(os.path.join(working_dir, fullname)):
            with open(os.path.join(working_dir, fullname), "r") as jcmfile:
                fobj[fullname] = jcmfile.read()
            break


GEOMETRY_PARAMS = {
    "sphere": ("radius",),
    "ellipsoid": ("radiusx", "radiusy", "radiusz"),
    "spheroid": ("radiusxy", "radiusz"),
    "cylinder": ("radius", "height"),
    "cone": ("radius", "height"),
    "torus": ("major_radius", "minor_radius"),
    "cube": ("length",),
    "rectangular_cuboid": ("lengthx", "lengthy", "lengthz"),
     "helix": ("radius_helix", "radius_wire", "number_turns", "pitch")
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
        fobj[param] = np.asarray(params[param])
    fobj.attrs["shape"] = str(shape)
    if meshfile is not None:
        mesh_data(fobj, meshfile, lunit)
        if isinstance(meshfile, str):
            with open(os.path.join(working_dir, meshfile), "rb") as fobj_mesh:
                fobj[os.path.split(meshfile)[1]] = np.void(fobj_mesh.read())

