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


def _check_name(fobj, name, exact=True):
    if (exact and fobj.name != name) or not fobj.name.startswith(name):
        raise ValueError(f"require name '{name}'")

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



